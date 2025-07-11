# Complete Workflow Example: From Sequence Design to Synthesis & Testing

This document demonstrates the complete workflow for designing sequences, selecting variants, and finding vendors for synthesis and testing.

## Workflow Overview

1. **Design sequences** by mutating a given sequence
2. **Select 96 variants** from the generated mutations
3. **Find synthesis vendors** for ordering the sequences
4. **Research testing options** for expression, function, and binding assays

---

## Step-by-Step Example

### Step 1: Design Sequences by Mutating a Given Sequence

**User Command:**
```
Mutate the following sequence and generate 96 variants:
ATGCGTACGTTAGCCTAGGCTAGCTAGGCTAGCTAGGCTAGCTAGGCTAGC
```

**System Response:**
- Generates 96 mutated variants
- Provides mutation statistics and visualization
- Creates downloadable CSV/FASTA file with all variants

### Step 2: Select 96 Variants

**User Command:**
```
Select all 96 variants for synthesis
```

**System Response:**
- Confirms selection of all 96 variants
- Prepares sequences for vendor research

### Step 3: Find Synthesis Vendors

**User Command:**
```
I want to order these 96 variants from a DNA synthesis vendor and test their expression, function, and binding
```

**System Response:**
```
DNA vendor research completed: Found 7 DNA synthesis vendors. Here are the options:

VENDORS FOUND:
1. Twist Bioscience
   - Services: DNA synthesis, Gene synthesis, Oligo pools
   - Pricing: $0.07-$0.09 per base pair
   - Turnaround: 5-10 business days
   - Max Length: Up to 3.2kb genes
   - Specialties: High-throughput, Gene libraries

2. GenScript
   - Services: Gene synthesis, DNA synthesis, Oligo synthesis
   - Pricing: $0.08-$0.12 per base pair
   - Turnaround: 7-14 business days
   - Max Length: Up to 10kb genes
   - Specialties: Custom genes, Codon optimization

3. Integrated DNA Technologies (IDT)
   - Services: Oligo synthesis, Gene synthesis, gBlocks
   - Pricing: $0.06-$0.10 per base pair
   - Turnaround: 2-5 business days
   - Max Length: Up to 2kb gBlocks
   - Specialties: High-quality oligos, Fast turnaround

4. Eurofins Genomics
   - Services: DNA synthesis, Gene synthesis, Oligo synthesis
   - Pricing: $0.08-$0.11 per base pair
   - Turnaround: 5-10 business days
   - Max Length: Up to 5kb genes
   - Specialties: High-throughput, Quality control

5. Synthego
   - Services: CRISPR gRNA synthesis, Gene synthesis, Oligo synthesis
   - Pricing: $0.05-$0.08 per base pair
   - Turnaround: 3-7 business days
   - Max Length: Up to 200bp oligos
   - Specialties: CRISPR tools, Gene editing

6. Thermo Fisher Scientific
   - Services: GeneArt gene synthesis, Oligo synthesis, Custom DNA
   - Pricing: $0.09-$0.13 per base pair
   - Turnaround: 7-14 business days
   - Max Length: Up to 15kb genes
   - Specialties: High-fidelity synthesis, Large genes

7. Biosearch Technologies
   - Services: Oligo synthesis, Custom DNA, Probe synthesis
   - Pricing: $0.07-$0.10 per base pair
   - Turnaround: 3-7 business days
   - Max Length: Up to 200bp oligos
   - Specialties: FISH probes, qPCR probes

RECOMMENDATIONS:
- For 96 variants, consider Twist Bioscience or Eurofins for high-throughput synthesis
- For fast turnaround, IDT provides 2-5 business days
- For large genes, GenScript or Thermo Fisher can handle up to 10-15kb
- Always request quotes from multiple vendors for best pricing
```

### Step 4: Research Testing Options

**User Command:**
```
What testing options are available for my sequences?
```

**System Response:**
```
DNA vendor research completed: Found 4 testing options. Here are the available assays:

TESTING OPTIONS:

1. Quality Control Testing
   - Services: Sanger sequencing, Mass spectrometry, HPLC analysis
   - Vendors: Eurofins Genomics, Genewiz, Macrogen
   - Pricing: $10-$50 per sample
   - Turnaround: 2-5 business days
   - Description: Verify sequence accuracy and purity

2. Validation Assays
   - Services: Functional testing, Expression analysis, Activity assays
   - Vendors: Charles River Labs, Eurofins Scientific, Covance
   - Pricing: $100-$500 per assay
   - Turnaround: 1-4 weeks
   - Description: Test biological function and activity

3. Stability Testing
   - Services: Temperature stability, pH stability, Long-term storage
   - Vendors: Stability Testing Labs, Eurofins, SGS
   - Pricing: $200-$1000 per study
   - Turnaround: 2-8 weeks
   - Description: Assess stability under various conditions

4. Purity Analysis
   - Services: HPLC, Mass spec, Gel electrophoresis
   - Vendors: Eurofins Genomics, Genewiz, Local labs
   - Pricing: $25-$100 per sample
   - Turnaround: 1-3 business days
   - Description: Determine purity and concentration

RECOMMENDATIONS:
- Start with Sanger sequencing for sequence verification
- Use HPLC or mass spec for purity analysis
- Consider functional assays for biological validation
- Always include quality control testing for any synthesis order
```

### Step 5: Vendor Comparison and Selection

**User Command:**
```
Compare Twist Bioscience and GenScript for my 96 variants
```

**System Response:**
```
VENDOR COMPARISON:

Twist Bioscience:
- Strengths: Large-scale synthesis, Gene libraries, High-throughput
- Limitations: Higher cost for small orders, Longer turnaround for custom work
- Best For: Gene libraries, Large-scale projects, Custom genes
- Contact: +1-650-938-6300, orders@twistbioscience.com
- Pricing: $0.07-$0.09 per base pair
- Turnaround: 5-10 business days (standard), 3-5 business days (rush)

GenScript:
- Strengths: Codon optimization, Custom genes, Mutagenesis
- Limitations: Higher pricing, Complex ordering process
- Best For: Custom gene synthesis, Codon optimization, Mutagenesis
- Contact: +1-732-885-9188, orders@genscript.com
- Pricing: $0.08-$0.12 per base pair
- Turnaround: 7-14 business days (standard), 5-7 business days (rush)

RECOMMENDATION:
For 96 variants, Twist Bioscience would be more cost-effective due to their high-throughput capabilities and lower per-base pricing. However, if you need codon optimization or custom modifications, GenScript might be better despite the higher cost.
```

---

## Natural Language Commands That Work

The system recognizes these commands and routes them to the appropriate tools:

### For Sequence Design:
- "Mutate this sequence"
- "Generate variants of this DNA"
- "Create mutations of the sequence"

### For Vendor Research:
- "I want to order the sequences from a DNA synthesis vendor"
- "Find DNA synthesis companies"
- "Research vendors for gene synthesis"
- "What vendors can synthesize my sequences?"

### For Testing Research:
- "What testing options are available for my sequences?"
- "Find assays to validate my sequences"
- "Research quality control testing for DNA"
- "What validation methods can I use?"

### For Specific Vendor Information:
- "Tell me about Twist Bioscience"
- "What are the details for GenScript?"
- "Compare Twist Bioscience and GenScript"

---

## Expected Output Format

The system provides structured responses with:
- Vendor information (name, services, pricing, turnaround)
- Testing options (assay types, vendors, pricing)
- Recommendations based on your requirements
- Contact information and ordering details
- Visualizations of vendor and testing options

This workflow enables researchers to seamlessly move from sequence design to experimental planning and vendor selection. 