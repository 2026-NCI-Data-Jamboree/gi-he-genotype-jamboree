
import copy
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor

SRC = '/home/lab/.claude-science/orgs/88dc4a54-20a0-46f3-ba12-f267bfbc3554/artifacts/proj_e3fbf9f86ff0/e54c2acb-397c-47f1-ae06-6d96d7dc75c5/v3a906c11_2026NCIDataJamboreeTemplateforReportOuts.pptx'
prs = Presentation(SRC)
s1,s2,s3,s4 = prs.slides[0],prs.slides[1],prs.slides[2],prs.slides[3]
NAVY = RGBColor(0x1F,0x3B,0x73); DARK = RGBColor(0x20,0x20,0x20); GREY=RGBColor(0x55,0x55,0x55)

def shp(slide,name):
    for s in slide.shapes:
        if s.name==name: return s
    raise KeyError(name)

def fill(tf, lines, default_size=16):
    """lines: list of dict(text, size, bold, color, level, space_after)"""
    tf.word_wrap = True
    tf.clear()
    for i,ln in enumerate(lines):
        p = tf.paragraphs[0] if i==0 else tf.add_paragraph()
        p.level = ln.get('level',0)
        p.space_after = Pt(ln.get('space_after',6))
        r = p.add_run(); r.text = ln['text']
        f = r.font
        f.name='Arial'; f.size=Pt(ln.get('size',default_size)); f.bold=ln.get('bold',False)
        f.color.rgb = ln.get('color',DARK)

def clone_slide(prs, src):
    new = prs.slides.add_slide(src.slide_layout)
    for shape in list(new.shapes):
        shape._element.getparent().remove(shape._element)
    for shape in src.shapes:
        new.shapes._spTree.append(copy.deepcopy(shape._element))
    return new

# ---------------- SLIDE 1 ----------------
tb = shp(s1,'TextBox 1')
tb.left,tb.top,tb.width,tb.height = Inches(2.95),Inches(1.62),Inches(9.9),Inches(1.25)
fill(tb.text_frame,[
 dict(text="Project Title: Reading molecular status from routine H&E across three NCI commons — "
           "a gastric and colorectal benchmark with pre-registration and two external cohorts", size=17, color=NAVY, bold=True),
 dict(text="Project Category: <<TEAM TO CONFIRM — e.g. multi-modal data integration / AI-ML>>", size=16, color=GREY),
])
tbl = shp(s1,'Table 4').table
team = [("<<Name>>","<<Institution>>","Project lead; design, analysis, interpretation"),
        ("<<Name>>","<<Institution>>","<<Role>>")]
for ri,(n,inst,role) in enumerate(team, start=1):
    for ci,val in enumerate((n,inst,role)):
        c = tbl.cell(ri,ci); c.text=val
        for p in c.text_frame.paragraphs:
            for r in p.runs: r.font.size=Pt(14); r.font.name='Arial'

# ---------------- SLIDE 2 ----------------
tb = shp(s2,'TextBox 1')
tb.left,tb.top,tb.width,tb.height = Inches(2.95),Inches(1.55),Inches(10.1),Inches(4.3)
fill(tb.text_frame,[
 dict(text="Specific problem(s) to solve:", size=19, bold=True, color=NAVY, space_after=3),
 dict(text="• MSI/mismatch-repair status and driver-mutation status decide immunotherapy and targeted therapy in gastric "
           "and colorectal cancer, but molecular testing is not universal and adds days to weeks — while the diagnostic "
           "H&E slide is already digitised and public in the Imaging Data Commons.", size=15, space_after=12),
 dict(text="Gaps/unmet needs to address:", size=19, bold=True, color=NAVY, space_after=3),
 dict(text="• Imaging (IDC), genomic (GDC) and proteomic (PDC) commons are not linked at the patient level, and published "
           "H&E-to-molecular models are rarely re-tested outside the commons they were trained in.", size=15, space_after=12),
 dict(text="Brief Description of overall approach (high-level):", size=19, bold=True, color=NAVY, space_after=3),
 dict(text="• Frozen pathology foundation-model features (UNI v1 and Phikon) with simple classifier heads over 1,396 GI "
           "whole-slide images from IDC; labels from GDC/cBioPortal; protein from PDC.", size=15, space_after=4),
 dict(text="• Pre-registered design: hospital-grouped 5-fold cross-validation, a locked held-out test, and two independent "
           "external cohorts (CPTAC-COAD and SurGen, UK).", size=15),
])
shp(s2,'TextBox 3')._element.getparent().remove(shp(s2,'TextBox 3')._element)

# ---------------- SLIDE 3 ----------------
t = shp(s3,'Table 2')
t.top, t.height, t.left, t.width = Inches(2.15), Inches(3.4), Inches(0.45), Inches(12.45)
tb3 = t.table
rows = [
 ("IDC: TCGA-COAD / READ / STAD (1,175 WSI) + CPTAC-COAD (221 WSI); 1,396 slides total",
  "Multi-modal: slides linked to GDC genomics by patient barcode; CPTAC arm also linked to PDC proteomics",
  "Tissue masking, tiling at 0.5 µm/px, frozen foundation-model embedding, mean-pool + linear / MLP / attention-MIL heads",
  "s5cmd, wsidicom, Phikon + UNI v1 encoders, PyTorch, scikit-learn; Modal T4 GPU workers",
  "Cloud"),
 ("GDC / cBioPortal: MSI (MANTIS), somatic mutation calls and TMB for 928 sequenced GI patients",
  "Genomic labels joined to the imaging cohort; patient-level, hospital-grouped splits",
  "Pre-registered analysis plan; 5-fold hospital-grouped CV; locked 20% test; single-look external reads logged as deviations",
  "pandas / scikit-learn; preregistration.md with spec hashes",
  "Cloud"),
 ("PDC000116 + PDC000109 (Prospective COAD Proteome, CPTAC2)",
  "Three-commons table: 94 patients with slide + genomic MSI call + MMR protein abundance",
  "MMR protein (MSH2 / MSH6 / MSH3 / EPCAM) vs image-derived MSI score; Mann-Whitney + Spearman",
  "PDC GraphQL API (geneSpectralCount, filesPerStudy); CPTAC protein assembly reports",
  "Cloud"),
 ("EBI BioStudies S-BIAD1285 (SurGen, UK colorectal; 128 CZI slides, 529 GB downloaded)",
  "Second independent external cohort; MMR-IHC and KRAS assay labels",
  "Frozen model applied once, no re-training; external MSI and MSS-restricted KRAS read-outs",
  "aria2c, pylibCZIrw, Modal volumes",
  "Cloud"),
 ("Public tools / models (no new data): HoVer-Net MoNuSAC, Phikon (HF), UNI v1 (gated HF)",
  "Image-derived nuclear morphometry fused with encoder scores",
  "Cell-type composition per slide; fusion vs encoder-only, paired bootstrap of the difference",
  "tiatoolbox 2.1.3, Hugging Face Hub",
  "Cloud"),
]
# clone rows as needed
tr_proto = tb3.rows[1]._tr
while len(tb3.rows) < len(rows)+1:
    tb3._tbl.append(copy.deepcopy(tr_proto))
for ri,row in enumerate(rows, start=1):
    for ci,val in enumerate(row):
        c = tb3.cell(ri,ci)
        c.text_frame.word_wrap=True
        c.text = val
        c.margin_top=Emu(18288); c.margin_bottom=Emu(18288)
        for p in c.text_frame.paragraphs:
            p.space_after=Pt(0)
            for r in p.runs: r.font.size=Pt(9.5); r.font.name='Arial'
for ci in range(5):
    c = tb3.cell(0,ci)
    for p in c.text_frame.paragraphs:
        for r in p.runs: r.font.size=Pt(11); r.font.name='Arial'
for w,ci in zip([2.75,2.55,2.95,2.9,0.9],range(5)):
    tb3.columns[ci].width = Inches(w)
nt = shp(s3,'TextBox 3'); nt.top=Inches(5.72); nt.left=Inches(0.45); nt.width=Inches(12.45)

# ---------------- SLIDE 4 : RESULTS 1 ----------------
ttl = shp(s4,'TextBox 8')
fill(ttl.text_frame,[dict(text="Results 1 — what we built, and what reproduced", size=30, bold=False, color=NAVY)])
ttl.left,ttl.top,ttl.width = Inches(0.5),Inches(0.45),Inches(12.4)
body = shp(s4,'TextBox 3')
body.left,body.top,body.width,body.height = Inches(0.5),Inches(1.35),Inches(12.4),Inches(5.6)
B=dict(size=14, space_after=9)
fill(body.text_frame,[
 dict(text="• Cohort-scale pipeline, run twice: all 1,396 GI slides from IDC embedded with two independent foundation "
           "models (1,395 with usable features), ~4,000 tissue tiles per slide — about $30 of cloud GPU for the whole project.", **B),
 dict(text="• MSI/dMMR reproduces the published benchmark with no task-specific training of the encoder: all-GI "
           "hospital-grouped CV 0.765 (Phikon) / 0.768 (UNI); locked TCGA test 0.760 / 0.826. Kather et al. 2019 reported "
           "0.77 (TCGA colorectal) and 0.81 (TCGA stomach) on the same source data.", **B),
 dict(text="• It holds on two cohorts outside TCGA: CPTAC-COAD, 140 slides — 0.701 / 0.794, matching the published "
           "CPTAC-COAD validation (Ding 2022, 0.73); SurGen UK, 126 slides with MMR immunohistochemistry labels — "
           "0.760 / 0.665, to our knowledge the first read-out of this model family on that newly released cohort.", **B),
 dict(text="• Three commons joined at the patient: 94 CPTAC patients carry slide + genomic MSI call + MMR protein. "
           "MSH2 and MSH6 protein are lower in MSI-H tumours (p = 0.010 / 0.031) and MSH6 abundance tracks the image MSI "
           "score (rho = -0.25, p = 0.014) — independent, non-imaging support that the model reads mismatch-repair biology.", **B),
 dict(text="• What that could mean in clinic (illustrative, single cohort): at the operating point giving 100% sensitivity "
           "on CPTAC-COAD, 12.9% of slides are ruled out of confirmatory MSI testing with no MSI-H case missed "
           "(35 MSI-H of 140 slides).", **B),
])
# remove the decorative left picture so full-width text fits
for s in list(s4.shapes):
    if s.shape_type==13: s._element.getparent().remove(s._element)

# ---------------- SLIDE 5 : '/home/lab/.claude-science/orgs/88dc4a54-20a0-46f3-ba12-f267bfbc3554/artifacts/proj_e3fbf9f86ff0/91921682-b6ea-4885-b7c3-19c40fed8735/ve32fcd3d_fig_forest_endpoints.png' FIGURE ----------------
s5 = clone_slide(prs, s4)
ttl5 = shp(s5,'TextBox 8')
fill(ttl5.text_frame,[dict(text="Results 2 — every endpoint we tested, on one axis", size=30, color=NAVY)])
b5 = shp(s5,'TextBox 3')
b5.left,b5.top,b5.width,b5.height = Inches(8.1),Inches(1.55),Inches(4.9),Inches(5.2)
fill(b5.text_frame,[
 dict(text="• One row per endpoint and evaluation; bars are 95% bootstrap confidence intervals; "
           "shaded rows are cohorts outside TCGA.", size=14, space_after=12),
 dict(text="• MSI clears chance in every evaluation and with both encoders, including both external cohorts — "
           "it is the endpoint that travels.", size=14, space_after=12),
 dict(text="• Colorectal TP53 and KRAS are learnable inside TCGA; outside it their intervals touch chance. "
           "Gastric CDH1 is the strongest in-TCGA signal but has no external cohort yet.", size=14, space_after=12),
 dict(text="• Encoder ranking is unresolved: UNI leads on CPTAC, Phikon on SurGen "
           "(paired difference -0.095, 95% CI -0.223 to 0.034).", size=14),
])
pic = s5.shapes.add_picture('/home/lab/.claude-science/orgs/88dc4a54-20a0-46f3-ba12-f267bfbc3554/artifacts/proj_e3fbf9f86ff0/91921682-b6ea-4885-b7c3-19c40fed8735/ve32fcd3d_fig_forest_endpoints.png', Inches(0.35), Inches(1.5), height=Inches(5.3))
pic.crop_top = 0.047
pic.width = int(Inches(5.3) * (1546/(1222*(1-0.047))))

# ---------------- SLIDE 6 : RESULTS 3 ----------------
s6 = clone_slide(prs, s4)
fill(shp(s6,'TextBox 8').text_frame,[dict(text="Results 3 — what the built-in controls added", size=30, color=NAVY)])
b6 = shp(s6,'TextBox 3')
b6.left,b6.top,b6.width,b6.height = Inches(0.5),Inches(1.3),Inches(12.4),Inches(5.7)
fill(b6.text_frame,[
 dict(text="• Tumour mutational burden from H&E is largely an MSI signal. Our MSI model, never trained on burden, reads "
           "TMB-high at 0.756 / 0.759 across ~720 GI slides — and 0.639 / 0.672 once MSI-H cases are removed. "
           "Decomposition: 142 of 171 TMB-high patients are MSI-H; 10 of the 29 MSS ones carry a POLE/POLD1 proofreading hotspot.", **B),
 dict(text="• Colorectal TP53 is learnable inside TCGA (CV 0.76, locked test 0.70) but is not confirmed on the first "
           "external test of this endpoint (CPTAC-COAD, 80 patients: 0.642 / 0.592). The signal is flat across truncating, "
           "hotspot-missense and other-missense variants, so it is not the p53-immunohistochemistry accumulation phenotype.", **B),
 dict(text="• Gastric CDH1 is more than Lauren histotype: 0.844 / 0.826 on 241 slides, and 0.896 / 0.867 within "
           "non-diffuse tumours alone, against 0.659 for histotype on its own. Based on 22 mutant cases, no external cohort yet.", **B),
 dict(text="• Nuclear morphometry (HoVer-Net, 1,337 slides): our per-slide cell compositions independently reproduce "
           "Zhang et al. 2025's colorectal-vs-gastric split (0.89 / 0.07 vs 0.65 / 0.27 epithelial / lymphocyte), but a "
           "six-feature nuclear summary added nothing to the encoders externally (0.701→0.722 and 0.794→0.807; both differences span zero).", **B),
 dict(text="• Pre-registered blinded reader arm on matched slide pairs: MSI 6 of 6 pairs called correctly; gastric CDH1 "
           "14 of 19 pairs, and 13 of 14 when all three readers agreed — readers recorded morphology scores and rationale before unblinding.", **B),
])

# ---------------- SLIDE 7 : FEEDBACK ----------------
s7 = clone_slide(prs, s4)
fill(shp(s7,'TextBox 8').text_frame,[dict(text="Feedback to the commons, pain points and open outputs", size=28, color=NAVY)])
b7 = shp(s7,'TextBox 3')
b7.left,b7.top,b7.width,b7.height = Inches(0.5),Inches(1.25),Inches(12.4),Inches(5.8)
S=dict(size=13, space_after=4)
fill(b7.text_frame,[
 dict(text="Gaps in the data we hit, and what each one cost us", size=16, bold=True, color=NAVY, space_after=5),
 dict(text="• PDC000116's quantitative matrix (7,418 genes) contains neither MLH1 nor PMS2, so the planned MLH1-loss check "
           "became an MSH2/MSH6 check; PDC000109 carries MLH1 but detects it in only 12 of 101 runs.", **S),
 dict(text="• No TCGA-matched stomach proteome exists — the protein leg is colorectal-only.", **S),
 dict(text="• TCGA-ESCA has 158 slides and zero MSI calls, so an oesophageal false-alarm control is impossible at any budget.", **S),
 dict(text="• 46 TCGA slides that carry an MSI call have no S3 URL in the manifests we could assemble, and were excluded.", space_after=11, **{k:v for k,v in S.items() if k!='space_after'}),
 dict(text="Interoperability and access pain points", size=16, bold=True, color=NAVY, space_after=5),
 dict(text="• PDC: quantDataMatrix returns 'Matrix data not found' for several studies; filesPerStudy returns empty records "
           "when queried by pdc_study_id but full metadata by study_id — undocumented, and a day of work to discover.", **S),
 dict(text="• Three commons, three patient-identifier conventions: the patient-level image–genomic–protein join had to be rebuilt by hand.", **S),
 dict(text="• No common slide format across sources — DICOM at 0.25–0.5 µm/px in IDC, Zeiss CZI at 0.11 µm/px in SurGen — "
           "each needing its own reader, resampling and tissue-masking path.", space_after=11, **{k:v for k,v in S.items() if k!='space_after'}),
 dict(text="Potential benefits to the broad community", size=16, bold=True, color=NAVY, space_after=5),
 dict(text="• Open release planned: pre-registration with spec hashes, pipeline code, per-slide scores and all result tables.", **S),
 dict(text="• Reusable assets: the 94-patient three-commons table, the cross-commons defect catalogue above, and an external-cohort "
           "options table for GI H&E work (which public cohorts carry which labels, and at what access cost).", **S),
])
prs.save('NCI_Jamboree_ReportOut_GI_HE_molecular.pptx')
print("saved", 'NCI_Jamboree_ReportOut_GI_HE_molecular.pptx')
