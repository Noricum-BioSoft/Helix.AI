import redis
import json
import asyncio
from typing import Dict, Any, Optional
import os
from contextlib import asynccontextmanager

# Redis connection
def get_redis_client():
    redis_url = os.getenv("REDIS_URL", "redis://localhost:6379")
    return redis.from_url(redis_url, decode_responses=True)

# Event bus utilities
class EventBus:
    def __init__(self):
        self.redis_client = get_redis_client()
    
    async def publish_event(self, channel: str, event_data: Dict[str, Any]):
        """Publish an event to a Redis channel"""
        try:
            self.redis_client.publish(channel, json.dumps(event_data))
        except Exception as e:
            print(f"Error publishing event: {e}")
    
    async def subscribe_to_events(self, channel: str, callback):
        """Subscribe to events from a Redis channel"""
        pubsub = self.redis_client.pubsub()
        pubsub.subscribe(channel)
        
        try:
            for message in pubsub.listen():
                if message["type"] == "message":
                    data = json.loads(message["data"])
                    await callback(data)
        except Exception as e:
            print(f"Error in event subscription: {e}")
        finally:
            pubsub.close()

# Database utilities
@asynccontextmanager
async def get_db_connection():
    """Get database connection context manager"""
    import asyncpg
    
    db_url = os.getenv("POSTGRES_URL", "postgresql://databloom:databloom123@localhost:5432/databloom")
    conn = await asyncpg.connect(db_url)
    try:
        yield conn
    finally:
        await conn.close()

# Service discovery
class ServiceRegistry:
    def __init__(self):
        self.redis_client = get_redis_client()
    
    async def register_service(self, service_name: str, service_url: str):
        """Register a service with the registry"""
        self.redis_client.hset("services", service_name, service_url)
    
    async def get_service_url(self, service_name: str) -> Optional[str]:
        """Get service URL from registry"""
        return self.redis_client.hget("services", service_name)
    
    async def list_services(self) -> Dict[str, str]:
        """List all registered services"""
        return self.redis_client.hgetall("services")

# Circuit breaker pattern
class CircuitBreaker:
    def __init__(self, failure_threshold: int = 5, timeout: int = 60):
        self.failure_threshold = failure_threshold
        self.timeout = timeout
        self.failure_count = 0
        self.last_failure_time = None
        self.state = "CLOSED"  # CLOSED, OPEN, HALF_OPEN
    
    async def call(self, func, *args, **kwargs):
        if self.state == "OPEN":
            if self.last_failure_time and (asyncio.get_event_loop().time() - self.last_failure_time) > self.timeout:
                self.state = "HALF_OPEN"
            else:
                raise Exception("Circuit breaker is OPEN")
        
        try:
            result = await func(*args, **kwargs)
            if self.state == "HALF_OPEN":
                self.state = "CLOSED"
                self.failure_count = 0
            return result
        except Exception as e:
            self.failure_count += 1
            self.last_failure_time = asyncio.get_event_loop().time()
            
            if self.failure_count >= self.failure_threshold:
                self.state = "OPEN"
            
            raise e

# Retry mechanism
async def retry_with_backoff(func, max_retries: int = 3, base_delay: float = 1.0):
    """Retry function with exponential backoff"""
    for attempt in range(max_retries):
        try:
            return await func()
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            
            delay = base_delay * (2 ** attempt)
            await asyncio.sleep(delay)

# Logging utilities
import logging

def setup_logging(service_name: str):
    """Setup logging for a service"""
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s - {service_name} - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(f'{service_name}.log')
        ]
    )
    return logging.getLogger(service_name)

# Metrics utilities
from prometheus_client import Counter, Histogram, Gauge

# Define metrics
request_counter = Counter('http_requests_total', 'Total HTTP requests', ['method', 'endpoint', 'status'])
request_duration = Histogram('http_request_duration_seconds', 'HTTP request duration')
active_workflows = Gauge('active_workflows', 'Number of active workflows')
active_sessions = Gauge('active_sessions', 'Number of active sessions')

def track_request_metrics(method: str, endpoint: str, status: int, duration: float):
    """Track request metrics"""
    request_counter.labels(method=method, endpoint=endpoint, status=status).inc()
    request_duration.observe(duration) 