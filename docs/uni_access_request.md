# Getting UNI access (MahmoodLab) — fill-in sheet

## 0. Hugging Face account email (do this FIRST)
The gate form requires that your **primary** HF email is your institutional address.
Settings -> Account -> Emails: add your hospital/university email, verify it, set it as primary.
@gmail / @hotmail / @qq as primary is the single most common cause of denial.

## 1. Accept the terms
- https://huggingface.co/MahmoodLab/UNI      <- ViT-L/16, 1024-d  (the one we will run)
- https://huggingface.co/MahmoodLab/UNI2-h   <- ViT-H/14, 1536-d  (accept too; free, 30 s)

Form fields — copy and edit the bracketed parts:

**Full name** (no abbreviations, no initials):
    [Firstname Lastname]

**Affiliation** (spelled out in full, no acronyms — "Department of Medicine,
Memorial ... Hospital", not "MSKCC"):
    [Department, Institution, City, Country]

**Type of affiliation:**
    Academic / non-profit research

**Official institutional email** (must equal your HF primary email):
    [you@institution.edu]

**Intended research use:**
    Non-commercial academic research. I am an internal medicine resident taking part in
    the NCI Jamboree. I will use UNI as a frozen tile-level feature extractor on public
    H&E whole-slide images from the NCI Imaging Data Commons (TCGA colorectal, gastric
    and oesophageal cohorts, plus CPTAC-COAD as an external validation set), and train
    attention-based multiple-instance learning classifiers on the extracted features to
    predict microsatellite instability status and KRAS mutation status. Results will be
    presented in an academic talk and, if publishable, in a peer-reviewed paper with
    proper attribution to the UNI paper. I will not redistribute the weights and will
    not use the model or anything derived from it commercially.

**Checkbox:** I agree not to distribute the model; other users in my organisation must
register individually.  -> tick it (it is a condition of the licence)

## 2. Create a read token
Settings -> Access Tokens -> Create new token

Simplest: token type **Read** (classic). Works on gated repos you have been granted.

If you pick **Fine-grained** instead, you MUST tick:
    [x] Read access to contents of all public gated repos you can access
A fine-grained token without that box returns 401 on UNI even after approval — this is
the usual "I was approved but it still says restricted" problem.

Copy the token (hf_...) immediately; HF shows it once.

## 3. Hand it to this workspace
Left sidebar -> Customize -> Credentials -> Add -> Hugging Face, paste the token.
Do not paste the token into the chat.

## Licence note
UNI is CC-BY-NC-ND 4.0: non-commercial academic use only, and the restriction extends to
models trained on UNI outputs. An academic talk and paper are fine; anything commercial
needs prior approval from the authors.
