# Safety and Privacy Policy

## Overview

This document defines safety, privacy, and dual-use policies for Helix.AI. These policies ensure responsible handling of sensitive data and prevent misuse of the system for harmful purposes.

---

## Human Data and Privacy

### Principle

**Treat all human genomic and clinical data as sensitive and personally identifying.**

Genomic data is unique to individuals and can be used to identify them even when traditional identifiers (names, addresses) are removed.

---

### Data Handling Requirements

#### 1. No Raw Identifiers in Outputs

**Requirement:** Never include personally identifying information (PII) in analysis outputs or summaries.

**PII includes:**
- Patient names
- Medical record numbers
- Dates of birth
- Addresses
- Phone numbers
- Email addresses
- Social security numbers
- Full ZIP codes (use first 3 digits only)

**Examples:**

❌ **Bad:**
```json
{
  "sample": "Patient John Smith, DOB 1980-05-15",
  "diagnosis": "Breast cancer",
  "tumor_mutations": [...]
}
```

✅ **Good:**
```json
{
  "sample": "Sample_001",
  "diagnosis": "Breast cancer",
  "tumor_mutations": [...]
}
```

#### 2. Aggregate, Don't Individualize

**Requirement:** Report aggregate statistics rather than individual-level data when possible.

**Examples:**

❌ **Bad:**
```
Patient A: 15 mutations
Patient B: 23 mutations
Patient C: 8 mutations
```

✅ **Good:**
```
Cohort summary:
- Mean mutations per sample: 15.3 (SD: 7.5)
- Range: 8-23 mutations
- n=3 samples
```

#### 3. Redact Sensitive Metadata

**Requirement:** Remove or obfuscate metadata that could identify individuals.

**Redaction rules:**

| Metadata Type | Action |
|---|---|
| Patient ID | Replace with anonymized ID (Sample_001) |
| Age | Report age range (40-45) or just >18/<18 |
| Location | Report only country or state, not city/ZIP |
| Dates | Report relative dates ("Day 0", "Week 4") |
| Rare variants | Group very rare variants to prevent identification |

#### 4. Access Control and Audit Logs

**Requirement:** Maintain access logs for all operations on human data.

**Logs must include:**
- User ID
- Timestamp
- Operation performed
- Data accessed
- Results returned

#### 5. Data Retention and Deletion

**Requirement:** Respect data retention policies and deletion requests.

**Policies:**
- Session data: Delete after 7 days of inactivity
- Uploaded human data: Delete within 30 days or upon user request
- Analysis results: Anonymize or delete after 90 days
- Audit logs: Retain for 1 year minimum

---

### Compliance Requirements

#### HIPAA (US Health Insurance Portability and Accountability Act)

**If handling US clinical data:**
- [ ] Implement HIPAA-compliant encryption (at rest and in transit)
- [ ] Maintain audit logs of all data access
- [ ] Provide data deletion upon request
- [ ] Sign Business Associate Agreements (BAAs) with healthcare partners
- [ ] Conduct regular HIPAA compliance audits

#### GDPR (EU General Data Protection Regulation)

**If handling EU resident data:**
- [ ] Obtain explicit consent for data processing
- [ ] Provide data portability (export user data on request)
- [ ] Implement "right to be forgotten" (delete data on request)
- [ ] Conduct Data Protection Impact Assessments (DPIAs)
- [ ] Appoint Data Protection Officer (DPO) if required

#### Other Regulations

- **Canada (PIPEDA)**: Personal Information Protection and Electronic Documents Act
- **Australia (Privacy Act 1988)**: Australian Privacy Principles
- **ELSI (Ethical, Legal, and Social Implications)**: Follow best practices for genomic research

---

## Dual-Use and Biosafety

### Principle

**Refuse requests that could enable harmful biological misuse while providing safe alternatives.**

---

### Categories of Harmful Requests

#### 1. Pathogen Enhancement

**Refuse:**
- Optimizing virulence factors
- Increasing transmissibility of pathogens
- Evading immune detection
- Antibiotic resistance engineering for harmful purposes

**Example refusal:**
```
I cannot help with optimizing pathogen virulence as this poses biosafety risks.

Instead, I can provide information on:
- Biosafety best practices for pathogen research
- Ethical guidelines for gain-of-function research
- Regulatory frameworks (NIH, WHO guidelines)
```

#### 2. Bioweapon Development

**Refuse:**
- Weaponization of biological agents
- Large-scale production methods for pathogens
- Delivery mechanisms for biological agents
- Evading detection systems

**Example refusal:**
```
I cannot assist with bioweapon development as this violates international law (Biological Weapons Convention).

Instead, I can provide information on:
- Biosecurity protocols
- Detection and surveillance systems
- International biosafety frameworks
```

#### 3. Unethical Human Genetic Engineering

**Refuse:**
- Heritable genome editing (germline modifications)
- Enhancement for non-medical purposes
- Genetic discrimination tools
- Eugenics-related analyses

**Example refusal:**
```
I cannot assist with heritable human genetic modifications as this raises serious ethical concerns.

Instead, I can provide information on:
- Ethical frameworks for human gene therapy
- Regulatory guidelines (FDA, EMA)
- Somatic cell therapy (non-heritable)
```

#### 4. Environmental Harm

**Refuse:**
- Engineering invasive species
- Disrupting ecosystems
- Bypassing ecological safeguards
- Creating organisms that cannot be contained

---

### Legitimate Research vs. Harmful Misuse

**How to distinguish:**

| Legitimate | Harmful |
|---|---|
| Understanding pathogen biology for therapeutic development | Increasing pathogen virulence without therapeutic justification |
| Studying antibiotic resistance to develop new drugs | Creating antibiotic-resistant pathogens for release |
| Modeling disease transmission for public health | Optimizing disease transmission for malicious purposes |
| Germline research in animal models | Unregulated human germline editing |
| Vaccine development | Evading vaccine protection |

**When in doubt:**
1. Check if research has clear therapeutic or public health benefit
2. Verify that research is conducted under appropriate oversight (IRB, IBC)
3. Ensure research follows published ethical guidelines
4. Consult biosafety experts if uncertain

---

### Safe Alternatives

When refusing a harmful request, always provide safe alternatives:

**Template:**
```
I cannot assist with [harmful request] as it poses [biosafety/ethical] risks.

Instead, I can help with:
1. [Safe alternative 1]
2. [Safe alternative 2]
3. [Resource for ethical guidance]
```

**Example:**
```
User: "How can I make this virus more transm