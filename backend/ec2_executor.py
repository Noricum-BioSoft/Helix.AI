"""
EC2 Executor for Tool Generator Agent

Manages EC2 instances with pre-installed bioinformatics tools and executes
generated code on these instances via SSH.
"""

import os
import boto3
import paramiko
import tempfile
import logging
import time
from typing import Dict, Any, Optional
from pathlib import Path
import subprocess
from botocore.exceptions import ClientError

logger = logging.getLogger(__name__)


class EC2Executor:
    """Manages EC2 instances for executing bioinformatics tool code."""
    
    def __init__(self):
        self.region = os.getenv('AWS_REGION', 'us-east-1')
        self.ec2_client = boto3.client('ec2', region_name=self.region)
        self.instance_id = os.getenv('HELIX_EC2_INSTANCE_ID')  # Pre-existing instance ID
        self.key_name = os.getenv('HELIX_EC2_KEY_NAME')  # EC2 key pair name
        self.key_file_path = os.getenv('HELIX_EC2_KEY_FILE')  # Path to private key file
        self.instance_type = os.getenv('HELIX_EC2_INSTANCE_TYPE', 't3.medium')
        self.ami_id = os.getenv('HELIX_EC2_AMI_ID')  # Custom AMI with bioinformatics tools
        self.security_group_id = os.getenv('HELIX_EC2_SECURITY_GROUP_ID')
        self.subnet_id = os.getenv('HELIX_EC2_SUBNET_ID')
        
        # Cache for instance details
        self._instance_cache: Optional[Dict[str, Any]] = None
        self._cache_time: float = 0
        self._cache_ttl = 60  # Cache for 60 seconds
    
    def get_or_create_instance(self) -> Optional[str]:
        """
        Get existing instance ID or create a new one.
        Returns instance ID or None if creation fails.
        """
        # If instance ID is provided via env var, use it
        if self.instance_id:
            logger.info(f"Using pre-configured EC2 instance: {self.instance_id}")
            instance_status = self._check_instance_status(self.instance_id)
            
            if instance_status == "running":
                return self.instance_id
            elif instance_status == "not_found":
                # Instance doesn't exist - check if we should auto-create
                # Accept 'true', '1', 'yes', 'on', or any non-empty value (except 'false')
                auto_create_val = os.getenv('HELIX_EC2_AUTO_CREATE', 'false').lower().strip()
                auto_create = auto_create_val in ('true', '1', 'yes', 'on', 'auto-create') or (auto_create_val and auto_create_val != 'false')
                if auto_create:
                    logger.warning(f"Instance {self.instance_id} does not exist. Auto-creating new instance...")
                    # Clear the invalid instance ID and create a new one
                    self.instance_id = None
                    return self._create_instance()
                else:
                    logger.error(
                        f"Instance {self.instance_id} does not exist. "
                        f"Set HELIX_EC2_AUTO_CREATE=true to auto-create a new instance, "
                        f"or update HELIX_EC2_INSTANCE_ID with a valid instance ID."
                    )
                    return None
            elif instance_status in ["stopped", "stopping"]:
                logger.warning(f"Instance {self.instance_id} is {instance_status}. Attempting to start...")
                if self._start_instance(self.instance_id):
                    return self.instance_id
                else:
                    logger.error(f"Failed to start instance {self.instance_id}")
                    return None
            else:
                logger.warning(f"Instance {self.instance_id} is in state: {instance_status}")
                return None
        
        # Try to find existing running instance with tag
        existing_instance = self._find_existing_instance()
        if existing_instance:
            logger.info(f"Found existing running instance: {existing_instance}")
            return existing_instance
        
        # Create new instance if none exists
        logger.info("No existing instance found, creating new EC2 instance...")
        return self._create_instance()
    
    def _find_existing_instance(self) -> Optional[str]:
        """Find existing running EC2 instance tagged for Helix.AI."""
        try:
            response = self.ec2_client.describe_instances(
                Filters=[
                    {'Name': 'tag:Project', 'Values': ['Helix.AI']},
                    {'Name': 'tag:Purpose', 'Values': ['BioinformaticsTools']},
                    {'Name': 'instance-state-name', 'Values': ['running', 'pending']}
                ]
            )
            
            for reservation in response['Reservations']:
                for instance in reservation['Instances']:
                    instance_id = instance['InstanceId']
                    state = instance['State']['Name']
                    if state == 'running':
                        return instance_id
                    elif state == 'pending':
                        logger.info(f"Instance {instance_id} is pending, waiting...")
                        if self._wait_for_instance(instance_id):
                            return instance_id
            return None
        except Exception as e:
            logger.error(f"Error finding existing instance: {e}")
            return None
    
    def _create_instance(self) -> Optional[str]:
        """Create a new EC2 instance with bioinformatics tools."""
        try:
            # Validate required configuration
            if not self.key_name:
                logger.error(
                    "Cannot create EC2 instance: HELIX_EC2_KEY_NAME is not set. "
                    "Set HELIX_EC2_KEY_NAME to your EC2 key pair name."
                )
                return None
            
            # Use custom AMI if provided, otherwise use Amazon Linux 2
            image_id = self.ami_id or self._get_default_ami()
            
            # Build launch configuration
            launch_config = {
                'ImageId': image_id,
                'MinCount': 1,
                'MaxCount': 1,
                'InstanceType': self.instance_type,
                'KeyName': self.key_name,  # Required for SSH access
                'TagSpecifications': [
                    {
                        'ResourceType': 'instance',
                        'Tags': [
                            {'Key': 'Project', 'Value': 'Helix.AI'},
                            {'Key': 'Purpose', 'Value': 'BioinformaticsTools'},
                            {'Key': 'Name', 'Value': 'Helix.AI-Bioinformatics-Tools'}
                        ]
                    }
                ]
            }
            
            # Add IAM instance profile for S3 access (if available)
            # Try to find or create an IAM role with S3 access
            iam_role_name = os.getenv('HELIX_EC2_IAM_ROLE_NAME')
            if iam_role_name:
                # Get instance profile ARN for the role
                try:
                    iam_client = boto3.client('iam', region_name=self.region)
                    response = iam_client.get_instance_profile(InstanceProfileName=iam_role_name)
                    launch_config['IamInstanceProfile'] = {'Arn': response['InstanceProfile']['Arn']}
                    logger.info(f"Using IAM role: {iam_role_name}")
                except Exception as e:
                    logger.warning(f"Could not use IAM role {iam_role_name}: {e}. Instance will need credentials configured manually.")
            
            # Add user data (base64 encoded)
            user_data = self._get_user_data_script()
            if user_data:
                launch_config['UserData'] = user_data
            
            if self.security_group_id:
                launch_config['SecurityGroupIds'] = [self.security_group_id]
            elif self.subnet_id:
                # Create security group if needed
                sg_id = self._get_or_create_security_group()
                if sg_id:
                    launch_config['SecurityGroupIds'] = [sg_id]
            
            if self.subnet_id:
                launch_config['SubnetId'] = self.subnet_id
            
            # Launch instance
            response = self.ec2_client.run_instances(**launch_config)
            instance_id = response['Instances'][0]['InstanceId']
            
            logger.info(f"Created EC2 instance: {instance_id}")
            
            # Wait for instance to be running
            if self._wait_for_instance(instance_id):
                return instance_id
            else:
                logger.error(f"Instance {instance_id} failed to start")
                return None
                
        except Exception as e:
            logger.error(f"Error creating EC2 instance: {e}", exc_info=True)
            return None
    
    def _get_default_ami(self) -> str:
        """Get default Amazon Linux 2 AMI for the region."""
        try:
            # Try SSM Parameter Store first (most reliable)
            try:
                ssm_client = boto3.client('ssm', region_name=self.region)
                parameter_name = '/aws/service/ami-amazon-linux-latest/amzn2-ami-hvm-x86_64-gp2'
                response = ssm_client.get_parameter(Name=parameter_name)
                ami_id = response['Parameter']['Value']
                logger.info(f"Found latest Amazon Linux 2 AMI via SSM: {ami_id}")
                return ami_id
            except Exception as ssm_error:
                logger.debug(f"SSM parameter lookup failed: {ssm_error}, trying describe_images...")
            
            # Fallback: Query EC2 and sort by creation date
            response = self.ec2_client.describe_images(
                Owners=['amazon'],
                Filters=[
                    {'Name': 'name', 'Values': ['amzn2-ami-hvm-*-x86_64-gp2']},
                    {'Name': 'state', 'Values': ['available']}
                ]
            )
            if response['Images']:
                # Sort by CreationDate (most recent first)
                images = sorted(
                    response['Images'],
                    key=lambda x: x.get('CreationDate', ''),
                    reverse=True
                )
                ami_id = images[0]['ImageId']
                logger.info(f"Found latest Amazon Linux 2 AMI via describe_images: {ami_id}")
                return ami_id
        except Exception as e:
            logger.warning(f"Could not find default AMI: {e}")
        
        # Last resort: Use SSM with explicit region
        try:
            ssm_client = boto3.client('ssm', region_name=self.region)
            parameter_name = '/aws/service/ami-amazon-linux-latest/amzn2-ami-hvm-x86_64-gp2'
            response = ssm_client.get_parameter(Name=parameter_name)
            return response['Parameter']['Value']
        except Exception:
            pass
        
        # Final fallback: Hardcoded AMI IDs (may become outdated)
        logger.warning(f"Using hardcoded fallback AMI - may be outdated. Region: {self.region}")
        fallback_amis = {
            'us-east-1': 'ami-0c55b159cbfafe1f0',
            'us-west-2': 'ami-0d1cd67c26f5fca19',
            'eu-west-1': 'ami-0c94855ba95b798c7'
        }
        return fallback_amis.get(self.region, 'ami-0c55b159cbfafe1f0')
    
    def _get_user_data_script(self) -> str:
        """Get user data script to install bioinformatics tools on instance startup.
        
        Note: boto3's run_instances will automatically base64 encode this.
        """
        return """#!/bin/bash
# Install bioinformatics tools on EC2 instance startup

# Update system
yum update -y

# Install basic dependencies
yum install -y wget curl git gcc gcc-c++ make

# Install Miniconda (using version compatible with GLIBC 2.26)
cd /tmp
# Use Miniconda3-py38_4.9.2 which works with GLIBC 2.26
wget -q https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh
chmod +x Miniconda3-py38_4.9.2-Linux-x86_64.sh
bash Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p /opt/conda
rm -f Miniconda3-py38_4.9.2-Linux-x86_64.sh

# Add conda to PATH
export PATH="/opt/conda/bin:$PATH"

# Initialize conda
/opt/conda/bin/conda init bash

# Add bioconda channels
/opt/conda/bin/conda config --add channels bioconda
/opt/conda/bin/conda config --add channels conda-forge
/opt/conda/bin/conda config --set channel_priority strict

# Install common bioinformatics tools and Python packages (this takes 10-15 minutes)
/opt/conda/bin/conda install -y bbtools samtools bcftools bwa bowtie2 fastqc biopython pandas numpy boto3 botocore

# Create symlinks for easy access
ln -sf /opt/conda/bin/bbmerge.sh /usr/local/bin/bbmerge.sh || true
ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools || true
ln -sf /opt/conda/bin/bcftools /usr/local/bin/bcftools || true
ln -sf /opt/conda/bin/bwa /usr/local/bin/bwa || true
ln -sf /opt/conda/bin/bowtie2 /usr/local/bin/bowtie2 || true
ln -sf /opt/conda/bin/fastqc /usr/local/bin/fastqc || true

# Create working directory
mkdir -p /opt/helix-tools
chmod 777 /opt/helix-tools

# Signal completion
echo "Bioinformatics tools installation completed" > /opt/helix-tools/install.log
"""
    
    def _get_or_create_security_group(self) -> Optional[str]:
        """Get or create security group for EC2 instances."""
        if self.security_group_id:
            return self.security_group_id
        
        try:
            # Try to find existing security group
            response = self.ec2_client.describe_security_groups(
                Filters=[
                    {'Name': 'tag:Project', 'Values': ['Helix.AI']},
                    {'Name': 'tag:Purpose', 'Values': ['BioinformaticsTools']}
                ]
            )
            
            if response['SecurityGroups']:
                return response['SecurityGroups'][0]['GroupId']
            
            # Create new security group
            vpc_id = self._get_default_vpc()
            if not vpc_id:
                logger.error("Could not find default VPC")
                return None
            
            response = self.ec2_client.create_security_group(
                GroupName='helix-bioinformatics-tools',
                Description='Security group for Helix.AI bioinformatics tools EC2 instances',
                VpcId=vpc_id
            )
            sg_id = response['GroupId']
            
            # Add SSH access (restrict to your IP in production)
            self.ec2_client.authorize_security_group_ingress(
                GroupId=sg_id,
                IpPermissions=[
                    {
                        'IpProtocol': 'tcp',
                        'FromPort': 22,
                        'ToPort': 22,
                        'IpRanges': [{'CidrIp': '0.0.0.0/0'}]  # Restrict in production!
                    }
                ]
            )
            
            # Tag the security group
            self.ec2_client.create_tags(
                Resources=[sg_id],
                Tags=[
                    {'Key': 'Project', 'Value': 'Helix.AI'},
                    {'Key': 'Purpose', 'Value': 'BioinformaticsTools'}
                ]
            )
            
            return sg_id
            
        except Exception as e:
            logger.error(f"Error creating security group: {e}")
            return None
    
    def _get_default_vpc(self) -> Optional[str]:
        """Get default VPC ID."""
        try:
            response = self.ec2_client.describe_vpcs(
                Filters=[{'Name': 'isDefault', 'Values': ['true']}]
            )
            if response['Vpcs']:
                return response['Vpcs'][0]['VpcId']
        except Exception as e:
            logger.error(f"Error finding default VPC: {e}")
        return None
    
    def _check_instance_running(self, instance_id: str) -> bool:
        """Check if instance is running."""
        status = self._check_instance_status(instance_id)
        return status == "running"
    
    def _check_instance_status(self, instance_id: str) -> str:
        """
        Check instance status.
        Returns: 'running', 'stopped', 'pending', 'stopping', 'not_found', or 'error'
        """
        try:
            response = self.ec2_client.describe_instances(InstanceIds=[instance_id])
            if response['Reservations'] and response['Reservations'][0]['Instances']:
                state = response['Reservations'][0]['Instances'][0]['State']['Name']
                return state
            return "not_found"
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', '')
            if error_code == 'InvalidInstanceID.NotFound':
                return "not_found"
            logger.error(f"Error checking instance status: {e}")
            return "error"
        except Exception as e:
            logger.error(f"Error checking instance status: {e}")
            return "error"
    
    def _start_instance(self, instance_id: str) -> bool:
        """Start a stopped instance."""
        try:
            self.ec2_client.start_instances(InstanceIds=[instance_id])
            return self._wait_for_instance(instance_id)
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', '')
            if error_code == 'InvalidInstanceID.NotFound':
                logger.error(f"Instance {instance_id} does not exist")
            else:
                logger.error(f"Error starting instance: {e}")
            return False
        except Exception as e:
            logger.error(f"Error starting instance: {e}")
            return False
    
    def _wait_for_instance(self, instance_id: str, max_wait: int = 300) -> bool:
        """Wait for instance to be running and ready."""
        logger.info(f"Waiting for instance {instance_id} to be ready...")
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            try:
                response = self.ec2_client.describe_instances(InstanceIds=[instance_id])
                if response['Reservations']:
                    instance = response['Reservations'][0]['Instances'][0]
                    state = instance['State']['Name']
                    
                    if state == 'running':
                        # Check if instance is ready (status checks)
                        status = instance.get('StateReason', {}).get('Message', '')
                        if 'running' in status.lower() or not status:
                            logger.info(f"Instance {instance_id} is running")
                            # Wait a bit more for SSH to be ready
                            time.sleep(30)
                            return True
                    
                    logger.info(f"Instance state: {state}, waiting...")
                    time.sleep(10)
            except Exception as e:
                logger.warning(f"Error checking instance status: {e}")
                time.sleep(10)
        
        logger.error(f"Instance {instance_id} did not become ready within {max_wait} seconds")
        return False
    
    def get_instance_details(self, instance_id: str) -> Optional[Dict[str, Any]]:
        """Get instance details including public IP and connection info."""
        # Check cache
        if (self._instance_cache and 
            self._instance_cache.get('InstanceId') == instance_id and
            time.time() - self._cache_time < self._cache_ttl):
            return self._instance_cache
        
        try:
            response = self.ec2_client.describe_instances(InstanceIds=[instance_id])
            if response['Reservations']:
                instance = response['Reservations'][0]['Instances'][0]
                key_name = instance.get('KeyName', 'N/A')
                details = {
                    'InstanceId': instance_id,
                    'PublicIpAddress': instance.get('PublicIpAddress'),
                    'PrivateIpAddress': instance.get('PrivateIpAddress'),
                    'State': instance['State']['Name'],
                    'InstanceType': instance['InstanceType'],
                    'KeyName': key_name
                }
                
                # Log key pair info for debugging
                logger.info(f"Instance {instance_id} details: IP={details.get('PublicIpAddress')}, KeyName={key_name}, Configured key_name={self.key_name}")
                
                # Update cache
                self._instance_cache = details
                self._cache_time = time.time()
                
                return details
        except Exception as e:
            logger.error(f"Error getting instance details: {e}")
        
        return None
    
    def execute_code_on_instance(
        self,
        code: str,
        instance_id: str,
        timeout: int = 300
    ) -> Dict[str, Any]:
        """
        Execute Python code on EC2 instance via SSH.
        
        Args:
            code: Python code to execute
            instance_id: EC2 instance ID
            timeout: Execution timeout in seconds
            
        Returns:
            Dictionary with execution results
        """
        instance_details = self.get_instance_details(instance_id)
        if not instance_details:
            return {
                "status": "error",
                "error": f"Could not get details for instance {instance_id}"
            }
        
        public_ip = instance_details.get('PublicIpAddress')
        if not public_ip:
            return {
                "status": "error",
                "error": f"Instance {instance_id} does not have a public IP address"
            }
        
        # Get SSH key
        key_file = self._get_ssh_key_file()
        if not key_file:
            return {
                "status": "error",
                "error": "SSH key file not configured. Set HELIX_EC2_KEY_FILE environment variable."
            }
        
        # Expand user path (~) and resolve absolute path
        key_file = os.path.expanduser(key_file)
        key_file = os.path.abspath(key_file)
        
        if not os.path.exists(key_file):
            return {
                "status": "error",
                "error": f"SSH key file not found: {key_file}. Set HELIX_EC2_KEY_FILE environment variable."
            }
        
        # Check file permissions (should be 600)
        stat_info = os.stat(key_file)
        mode = oct(stat_info.st_mode)[-3:]
        if mode != '600':
            logger.warning(f"SSH key file {key_file} has permissions {mode}, should be 600. Attempting to fix...")
            try:
                os.chmod(key_file, 0o600)
                logger.info(f"Fixed permissions on {key_file}")
            except Exception as e:
                logger.warning(f"Could not fix permissions on {key_file}: {e}")
        
        # Verify key pair matches instance (if we have key_name info)
        instance_key_name = instance_details.get('KeyName')
        if instance_key_name and self.key_name and instance_key_name != self.key_name:
            logger.warning(
                f"Key pair mismatch: Instance {instance_id} uses key '{instance_key_name}', "
                f"but HELIX_EC2_KEY_NAME is set to '{self.key_name}'. "
                f"SSH authentication may fail."
            )
        elif instance_key_name:
            logger.info(f"Instance {instance_id} uses key pair: {instance_key_name}")
        
        try:
            # Create temporary Python file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
                f.write(code)
                local_script = f.name
            
            try:
                # Copy script to instance
                remote_script = f"/tmp/helix_exec_{int(time.time())}.py"
                self._scp_file(local_script, remote_script, public_ip, key_file)
                
                # Find Python executable on instance (try conda first, then system Python)
                python_cmd = self._find_python_on_instance(public_ip, key_file)
                logger.info(f"Using Python: {python_cmd}")
                
                # Ensure required Python packages are installed
                self._ensure_python_packages(public_ip, key_file, python_cmd)
                
                # Execute script on instance
                # EC2 instances automatically get AWS credentials from IAM role via instance metadata
                # Set AWS region environment variable for boto3
                # Include user site-packages in PYTHONPATH for pip-installed packages
                env_vars = (
                    f"AWS_DEFAULT_REGION={self.region} "
                    f"AWS_REGION={self.region} "
                    f"PYTHONPATH=/home/ec2-user/.local/lib/python3.8/site-packages:${{PYTHONPATH:-}}"
                )
                result = self._ssh_execute(
                    f"cd /opt/helix-tools && {env_vars} {python_cmd} {remote_script}",
                    public_ip,
                    key_file,
                    timeout
                )
                
                # Clean up remote script
                try:
                    self._ssh_execute(f"rm -f {remote_script}", public_ip, key_file, 10)
                except:
                    pass
                
                return result
                
            finally:
                # Clean up local script
                try:
                    os.unlink(local_script)
                except:
                    pass
                    
        except Exception as e:
            error_msg = str(e)
            error_type = type(e).__name__
            logger.error(f"Error executing code on EC2 instance: {error_type}: {error_msg}", exc_info=True)
            return {
                "status": "error",
                "error": f"{error_type}: {error_msg}",
                "error_type": error_type,
                "details": f"Failed to execute code on EC2 instance {instance_id}. Check logs for details."
            }
    
    def _get_ssh_key_file(self) -> Optional[str]:
        """Get path to SSH private key file."""
        if self.key_file_path:
            # Expand user path (~) and check if exists
            expanded = os.path.expanduser(self.key_file_path)
            if os.path.exists(expanded):
                return expanded
            # Return expanded path even if doesn't exist (will be checked later with better error)
            return expanded
        
        # Try common locations
        if self.key_name:
            common_paths = [
                f"~/.ssh/{self.key_name}.pem",
                f"~/.ssh/{self.key_name}",
                f"~/.ssh/id_rsa",
            ]
            for path in common_paths:
                expanded = os.path.expanduser(path)
                if os.path.exists(expanded):
                    return expanded
        
        return None
    
    def _scp_file(self, local_path: str, remote_path: str, host: str, key_file: str):
        """Copy file to EC2 instance via SCP."""
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        try:
            # Try to load private key (supports both RSA and ED25519)
            key = None
            key_errors = []
            
            # Try RSA first
            try:
                key = paramiko.RSAKey.from_private_key_file(key_file)
                logger.debug(f"Successfully loaded RSA key from {key_file}")
            except Exception as e:
                key_errors.append(f"RSA: {str(e)}")
                # Try Ed25519
                try:
                    key = paramiko.Ed25519Key.from_private_key_file(key_file)
                    logger.debug(f"Successfully loaded Ed25519 key from {key_file}")
                except Exception as e2:
                    key_errors.append(f"Ed25519: {str(e2)}")
                    # Try RSA with explicit password=None
                    try:
                        key = paramiko.RSAKey.from_private_key_file(key_file, password=None)
                        logger.debug(f"Successfully loaded RSA key (with password=None) from {key_file}")
                    except Exception as e3:
                        key_errors.append(f"RSA (password=None): {str(e3)}")
                        raise ValueError(
                            f"Failed to load SSH key from {key_file}. "
                            f"Tried RSA, Ed25519, and RSA with password=None. "
                            f"Errors: {'; '.join(key_errors)}"
                        )
            
            if not key:
                raise ValueError("Failed to load SSH key: key is None")
            
            # Connect
            logger.debug(f"Connecting to {host} as ec2-user with key from {key_file}")
            ssh.connect(
                host,
                username='ec2-user',
                pkey=key,
                timeout=30,
                look_for_keys=False,
                allow_agent=False
            )
            logger.debug(f"Successfully connected to {host}")
            
            # Use SFTP to copy file
            sftp = ssh.open_sftp()
            sftp.put(local_path, remote_path)
            sftp.close()
            logger.debug(f"Successfully copied {local_path} to {host}:{remote_path}")
            
        except paramiko.AuthenticationException as e:
            logger.error(f"SSH authentication failed for {host}: {e}")
            raise
        except paramiko.SSHException as e:
            logger.error(f"SSH error connecting to {host}: {e}")
            raise
        except Exception as e:
            logger.error(f"Error copying file to {host}: {e}", exc_info=True)
            raise
        finally:
            try:
                ssh.close()
            except:
                pass
    
    def _ssh_execute(self, command: str, host: str, key_file: str, timeout: int = 300) -> Dict[str, Any]:
        """Execute command on EC2 instance via SSH."""
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        try:
            # Try to load private key (supports both RSA and ED25519)
            key = None
            key_errors = []
            
            # Try RSA first
            try:
                key = paramiko.RSAKey.from_private_key_file(key_file)
                logger.debug(f"Successfully loaded RSA key from {key_file}")
            except Exception as e:
                key_errors.append(f"RSA: {str(e)}")
                # Try Ed25519
                try:
                    key = paramiko.Ed25519Key.from_private_key_file(key_file)
                    logger.debug(f"Successfully loaded Ed25519 key from {key_file}")
                except Exception as e2:
                    key_errors.append(f"Ed25519: {str(e2)}")
                    # Try RSA with explicit password=None
                    try:
                        key = paramiko.RSAKey.from_private_key_file(key_file, password=None)
                        logger.debug(f"Successfully loaded RSA key (with password=None) from {key_file}")
                    except Exception as e3:
                        key_errors.append(f"RSA (password=None): {str(e3)}")
                        raise ValueError(
                            f"Failed to load SSH key from {key_file}. "
                            f"Tried RSA, Ed25519, and RSA with password=None. "
                            f"Errors: {'; '.join(key_errors)}"
                        )
            
            if not key:
                raise ValueError("Failed to load SSH key: key is None")
            
            # Connect with retry
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    ssh.connect(
                        host,
                        username='ec2-user',
                        pkey=key,
                        timeout=30,
                        look_for_keys=False,
                        allow_agent=False
                    )
                    break
                except Exception as e:
                    if attempt < max_retries - 1:
                        logger.warning(f"SSH connection attempt {attempt + 1} failed, retrying...")
                        time.sleep(10)
                    else:
                        raise
            
            # Execute command
            stdin, stdout, stderr = ssh.exec_command(command, timeout=timeout)
            
            # Read output
            stdout_text = stdout.read().decode('utf-8')
            stderr_text = stderr.read().decode('utf-8')
            exit_status = stdout.channel.recv_exit_status()
            
            if exit_status == 0:
                return {
                    "status": "success",
                    "stdout": stdout_text,
                    "stderr": stderr_text,
                    "returncode": exit_status
                }
            else:
                # Command failed - include error message
                error_msg = stderr_text.strip() if stderr_text.strip() else f"Command failed with exit code {exit_status}"
                return {
                    "status": "error",
                    "error": error_msg,
                    "stdout": stdout_text,
                    "stderr": stderr_text,
                    "returncode": exit_status
                }
            
        except paramiko.SSHException as e:
            return {
                "status": "error",
                "error": f"SSH error: {str(e)}"
            }
        except Exception as e:
            return {
                "status": "error",
                "error": str(e)
            }
        finally:
            try:
                ssh.close()
            except:
                pass
    
    def _find_python_on_instance(self, host: str, key_file: str) -> str:
        """
        Find Python executable on EC2 instance.
        Prioritizes conda Python (3.8+) over system Python (3.7).
        """
        # Try conda Python first (preferred - Python 3.8+)
        conda_paths = [
            '/opt/conda/bin/python3',
            '/opt/conda/bin/python'
        ]
        
        for python_path in conda_paths:
            try:
                # Check if Python exists and get version
                result = self._ssh_execute(
                    f"test -f {python_path} && {python_path} --version 2>&1",
                    host,
                    key_file,
                    timeout=10
                )
                
                if result.get("status") == "success" and result.get("stdout", "").strip():
                    version_output = result.get("stdout", "").strip()
                    # Check if it's Python 3.8 or higher
                    import re
                    version_match = re.search(r'Python (\d+)\.(\d+)', version_output)
                    if version_match:
                        major, minor = int(version_match.group(1)), int(version_match.group(2))
                        if major == 3 and minor >= 8:
                            logger.info(f"Found conda Python {major}.{minor} at: {python_path}")
                            return python_path
                        else:
                            logger.warning(f"Conda Python at {python_path} is {major}.{minor}, need 3.8+")
            except Exception as e:
                logger.debug(f"Failed to check {python_path}: {e}")
                continue
        
        # Fallback to system Python (but warn if it's 3.7)
        system_paths = [
            '/usr/bin/python3',
            '/usr/bin/python',
            'python3',
            'python'
        ]
        
        for python_path in system_paths:
            try:
                # First find the path, then check version
                which_result = self._ssh_execute(
                    f"which {python_path} 2>/dev/null || command -v {python_path} 2>/dev/null",
                    host,
                    key_file,
                    timeout=10
                )
                
                if which_result.get("status") == "success" and which_result.get("stdout", "").strip():
                    actual_path = which_result.get("stdout", "").strip().split('\n')[0]
                    
                    # Now check version
                    version_result = self._ssh_execute(
                        f"{actual_path} --version 2>&1",
                        host,
                        key_file,
                        timeout=10
                    )
                    
                    if version_result.get("status") == "success":
                        version_output = version_result.get("stdout", "").strip()
                        import re
                        version_match = re.search(r'Python (\d+)\.(\d+)', version_output)
                        if version_match:
                            major, minor = int(version_match.group(1)), int(version_match.group(2))
                            if major == 3 and minor < 8:
                                logger.warning(f"Using Python {major}.{minor} at {actual_path} - boto3 may have deprecation warnings")
                            else:
                                logger.info(f"Found system Python {major}.{minor} at: {actual_path}")
                            return actual_path
            except Exception as e:
                logger.debug(f"Failed to check {python_path}: {e}")
                continue
        
        # Last resort: use python3 (should be available on most Linux systems)
        logger.warning("Could not verify Python path, using 'python3' as fallback")
        return "python3"
    
    def _ensure_python_packages(self, host: str, key_file: str, python_cmd: str):
        """
        Ensure required Python packages (boto3, etc.) are installed.
        Tries pip install if packages are missing.
        """
        required_packages = ['boto3', 'botocore']
        
        for package in required_packages:
            try:
                # Check if package is installed
                result = self._ssh_execute(
                    f"{python_cmd} -c 'import {package}' 2>&1",
                    host,
                    key_file,
                    timeout=10
                )
                
                if result.get("returncode") != 0:
                    logger.warning(f"Package {package} not found, attempting to install...")
                    
                    # Try pip install (works for both conda and system Python)
                    install_result = self._ssh_execute(
                        f"{python_cmd} -m pip install {package} --quiet 2>&1",
                        host,
                        key_file,
                        timeout=60
                    )
                    
                    if install_result.get("returncode") == 0:
                        logger.info(f"Successfully installed {package}")
                    else:
                        logger.warning(f"Failed to install {package}: {install_result.get('stderr', 'Unknown error')}")
                else:
                    logger.debug(f"Package {package} is already installed")
            except Exception as e:
                logger.warning(f"Could not check/install {package}: {e}")


# Global instance
_ec2_executor: Optional[EC2Executor] = None


def get_ec2_executor() -> EC2Executor:
    """Get or create EC2 executor instance."""
    global _ec2_executor
    if _ec2_executor is None:
        _ec2_executor = EC2Executor()
    return _ec2_executor

