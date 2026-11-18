
import copy
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from PIL import Image

prs = Presentation('/home/lab/.claude-science/orgs/88dc4a54-20a0-46f3-ba12-f267bfbc3554/artifacts/proj_e3fbf9f86ff0/f03e1a51-c528-49be-9061-48272e8fb320/v2af55bca_2026NCIDataJamboreeTemplateforReportOuts1.pptx')
NAVY=RGBColor(0x1F,0x3B,0x73); DARK=RGBColor(0x26,0x26,0x26); GREY=RGBColor(0x5A,0x5A,0x5A)
WHITE=RGBColor(0xFF,0xFF,0xFF); CARD=RGBColor(0xE8,0xEE,0xF6); ACC=RGBColor(0xC0,0x39,0x2B)
s1,s2,s3,s4 = prs.slides[0],prs.slides[1],prs.slides[2],prs.slides[3]
def shp(sl,n):
    for s in sl.shapes:
        if s.name==n: return s
    raise KeyError(n)
def rm(sh): sh._element.getparent().remove(sh._element)
def fill(tf,lines,wrap=True):
    tf.word_wrap=wrap; tf.clear()
    for i,ln in enumerate(lines):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph()
        p.space_after=Pt(ln.get('sa',6)); p.level=ln.get('level',0)
        if ln.get('align'): p.alignment=ln['align']
        if ln.get('lead'):
            r=p.add_run(); r.text=ln['lead']; r.font.name='Arial'; r.font.size=Pt(ln.get('size',16))
            r.font.bold=True; r.font.color.rgb=ln.get('leadcolor',NAVY)
        r=p.add_run(); r.text=ln['text']; f=r.font
        f.name='Arial'; f.size=Pt(ln.get('size',16)); f.bold=ln.get('bold',False); f.color.rgb=ln.get('color',DARK)
def clone(src):
    new=prs.slides.add_slide(src.slide_layout)
    for s in list(new.shapes): rm(s)
    for s in src.shapes: new.shapes._spTree.append(copy.deepcopy(s._element))
    return new
def content_slide(title, takeaway=None):
    sl=clone(s4)
    for s in list(sl.shapes):
        if s.shape_type==13: rm(s)
    t=shp(sl,'TextBox 8'); t.left,t.top,t.width,t.height=Inches(0.55),Inches(0.30),Inches(12.2),Inches(0.8)
    fill(t.text_frame,[dict(text=title,size=28,color=NAVY)])
    b=shp(sl,'TextBox 3'); b.left,b.top,b.width,b.height=Inches(0.55),Inches(1.08),Inches(12.2),Inches(0.6)
    if takeaway: fill(b.text_frame,[dict(text=takeaway,size=16,color=GREY)])
    else: rm(b)
    return sl
def put_image(sl,path,x,y,w,h):
    iw,ih=Image.open(path).size; sc=min(w/iw,h/ih)
    ww,hh=iw*sc,ih*sc
    return sl.shapes.add_picture(path, Inches(x+(w-ww)/2), Inches(y+(h-hh)/2), width=Inches(ww), height=Inches(hh))

# ---------- 1 Project Description ----------
t1=shp(s1,'TextBox 6'); t1.left,t1.width=Inches(3.9),Inches(6.6)
tb=shp(s1,'TextBox 1'); tb.left,tb.top,tb.width,tb.height=Inches(2.95),Inches(1.60),Inches(10.0),Inches(1.2)
fill(tb.text_frame,[
 dict(text="Reading the Genome Off the Slide — an interpretable pathology platform across NCI Cancer Data Commons",size=18,color=NAVY,bold=True,sa=4),
 dict(lead="Category:  ",text="data interoperability · cohort building · pipeline development · visualization",size=14,color=GREY)])
lbl=shp(s1,'TextBox 7'); lbl.left,lbl.top,lbl.width=Inches(2.95),Inches(2.95),Inches(9.0)
rm(shp(s1,'Table 4'))
au=s1.shapes.add_textbox(Inches(2.95),Inches(3.45),Inches(10.1),Inches(2.6))
fill(au.text_frame,[
 dict(text="Maria Alejandra Molina-Rodriguez, MD\u00b9   ·   Diksha Pasnoor, MBBS\u00b2   ·   Rabe Alhurani, MD\u00b9   ·   "
           "Furkan Haney\u00b9   ·   Meenal Gehlawat, MD\u00b9   ·   Ivonne Grande, MD\u00b9   ·   Dani Castillo, MD\u00b3",
      size=15,color=DARK,sa=10),
 dict(text="\u00b9 Desert Valley Hospital    \u00b2 Kamineni Academy of Medical Sciences and Research Centre, India    \u00b3 City of Hope",
      size=12,color=GREY,sa=14),
 dict(lead="Roles.  ",text="Project lead: M. A. Molina-Rodriguez.   Senior author, oncology oversight: D. Castillo.   "
           "Methodology: R. Alhurani.   Pipeline and training: F. Haney.   Cohort and validation: M. Gehlawat, I. Grande.   "
           "Documentation and figures: D. Pasnoor.", size=12,color=GREY,sa=0)])

# ---------- 2 Rationale ----------
t2=shp(s2,'TextBox 8'); t2.left,t2.top,t2.width=Inches(0.9),Inches(0.5),Inches(12.0)
bd=shp(s2,'TextBox 1'); bd.left,bd.top,bd.width,bd.height=Inches(3.1),Inches(2.05),Inches(9.9),Inches(4.0)
fill(bd.text_frame,[
 dict(lead="The problem.  ",text="MSI and mutation status decide GI therapy — but testing is uneven and slow.",size=20,sa=20),
 dict(lead="The gap.  ",text="Image, genomic and protein commons are not joined at the patient level.",size=20,sa=20),
 dict(lead="Our approach.  ",text="Frozen foundation models on 1,396 public slides. Pre-registered. Two external cohorts.",size=20,sa=6)])
rm(shp(s2,'TextBox 3'))

# ---------- 3 How it works ----------
sl=content_slide("How a slide becomes a number","The encoder is frozen — only the last, small step is fitted.")
put_image(sl,'fig_chain.png',0.5,2.6,12.3,1.7)
sl.notes_slide.notes_text_frame.text="No training of the image encoder. One logistic layer on top of frozen features."

# ---------- 4 Cohort in five numbers ----------
sl=content_slide("The cohort, in five numbers","Gastric, colorectal and rectal cancer. All public data.")
cards=[("1,396","whole-slide images\nImaging Data Commons"),("928","sequenced patients\nGenomic Data Commons"),
       ("3","NCI commons joined\nimage · genome · protein"),("2","independent external\ncohorts (US, UK)"),
       ("7","molecular endpoints\ntested")]
cw,gap=2.28,0.22; x0=(13.333-(len(cards)*cw+(len(cards)-1)*gap))/2
for i,(num,lab) in enumerate(cards):
    box=sl.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(x0+i*(cw+gap)),Inches(2.3),Inches(cw),Inches(2.5))
    box.fill.solid(); box.fill.fore_color.rgb=CARD; box.line.color.rgb=RGBColor(0xC3,0xD0,0xE4); box.line.width=Pt(1)
    box.shadow.inherit=False
    tf=box.text_frame; tf.word_wrap=True; tf.vertical_anchor=MSO_ANCHOR.MIDDLE
    fill(tf,[dict(text=num,size=40,bold=True,color=NAVY,align=PP_ALIGN.CENTER,sa=4)]+
            [dict(text=l,size=12,color=GREY,align=PP_ALIGN.CENTER,sa=0) for l in lab.split("\n")])

# ---------- 5 Data and tools (template table) ----------
rm(shp(s3,'TextBox 8'))
dt=shp(s3,'TextBox 4'); dt.left,dt.top,dt.width=Inches(0.55),Inches(0.30),Inches(12.2)
fill(dt.text_frame,[dict(text="Data and Tools Used",size=28,color=NAVY)])
tt=shp(s3,'Table 2'); tt.top,tt.left,tt.width,tt.height=Inches(1.30),Inches(0.5),Inches(12.3),Inches(4.4)
TB=tt.table
ROWS=[("IDC — TCGA-COAD / READ / STAD (1,175) + CPTAC-COAD (221)","Slides + GDC genomics + PDC proteomics, by patient",
       "Tile at 0.5 µm/px, frozen-encoder embedding, linear / attention-MIL head","Phikon, UNI v1, wsidicom, s5cmd, PyTorch","Cloud"),
      ("GDC / cBioPortal — MSI, mutations, TMB (928 patients)","Genomic labels joined to the slide cohort",
       "Hospital-grouped 5-fold CV, locked test, pre-registered","scikit-learn, preregistration.md","Cloud"),
      ("PDC000116 + PDC000109 — prospective COAD proteome","94 patients with slide + genome + protein",
       "MMR protein vs image score (Mann-Whitney, Spearman)","PDC GraphQL API","Cloud"),
      ("BioStudies S-BIAD1285 — SurGen UK colorectal (128 slides)","Second external cohort, MMR-IHC labels",
       "Frozen model applied once, no re-training","pylibCZIrw, aria2c, Modal","Cloud"),
      ("Public models — HoVer-Net MoNuSAC","Nuclear morphometry fused with encoder scores",
       "Cell-type composition per slide, paired bootstrap","tiatoolbox 2.1.3","Cloud")]
proto=TB.rows[1]._tr
while len(TB.rows) < len(ROWS)+1: TB._tbl.append(copy.deepcopy(proto))
for ri,rowv in enumerate(ROWS,start=1):
    for ci,v in enumerate(rowv):
        c=TB.cell(ri,ci); c.text_frame.word_wrap=True; c.text=v
        c.margin_top=Emu(27432); c.margin_bottom=Emu(27432)
        for p in c.text_frame.paragraphs:
            p.space_after=Pt(0)
            for r in p.runs: r.font.size=Pt(10.5); r.font.name='Arial'
for ci in range(5):
    for p in TB.cell(0,ci).text_frame.paragraphs:
        for r in p.runs: r.font.size=Pt(11); r.font.name='Arial'
for w_,ci in zip([3.1,2.7,3.2,2.5,0.8],range(5)): TB.columns[ci].width=Inches(w_)
nt=shp(s3,'TextBox 3'); nt.left,nt.top,nt.width=Inches(0.5),Inches(6.0),Inches(12.3)

# ---------- 6..11 results ----------
sl=content_slide("Does MSI detection hold up?","Lands in the published benchmark range — and stays there in two hospitals it never saw.")
put_image(sl,'fig_msi_scoreboard.png',0.5,1.75,12.3,5.2)
sl.notes_slide.notes_text_frame.text="Colorectal. Kather 2019: 0.77 TCGA, 0.84 external DACHS. Ding 2022 reported 0.73 on CPTAC-COAD."

sl=content_slide("A check that does not use the image at all","Mismatch-repair protein, measured by mass spectrometry, moves with the image score.")
put_image(sl,'fig_protein_mmr.png',2.2,1.75,8.9,5.1)
sl.notes_slide.notes_text_frame.text="94 CPTAC patients across three commons. MSH6 also correlates with the image MSI score: rho -0.25, p 0.014."

sl=content_slide("What this would mean in clinic","One in eight slides spared a confirmatory test — and no MSI-H case missed.")
put_image(sl,'fig_triage_waffle.png',1.2,1.8,10.9,5.0)
sl.notes_slide.notes_text_frame.text="Illustrative, single cohort, 35 positives. MSIntuit CE-IVD spares about one in two; we spare about one in eight."

sl=content_slide("Every endpoint we tested, on one axis","Shaded rows are cohorts outside TCGA.")
put_image(sl,'fig_all_endpoints.png',0.5,1.72,12.3,5.3)
sl.notes_slide.notes_text_frame.text="[DROPPABLE if short on time] Bars are 95% bootstrap intervals."

sl=content_slide("High mutational burden is mostly just MSI","Remove the MSI-H cases and the burden signal falls towards chance.")
put_image(sl,'fig_tmb.png',2.4,1.75,8.5,5.1)
sl.notes_slide.notes_text_frame.text="142 of 171 TMB-high patients are MSI-H; 10 of the 29 MSS ones carry a POLE/POLD1 hotspot."

sl=content_slide("What travels to another hospital — and what does not","MSI keeps its signal outside TCGA. TP53 and KRAS do not.")
put_image(sl,'fig_transfer.png',0.6,1.75,12.1,5.2)
sl.notes_slide.notes_text_frame.text="Gastric CDH1 is the strongest in-TCGA signal (0.84) but has no external cohort yet."

# ---------- 12 feedback ----------
sl=content_slide("Feedback to the commons","Everything below is dated, reproducible and in the written report.")
heads=[("Data gaps we hit",["PDC000116: no MLH1, no PMS2","No TCGA-matched stomach proteome","TCGA-ESCA: 158 slides, zero MSI calls","46 labelled slides with no file URL"]),
       ("Interoperability pain points",["Three commons, three patient-ID schemes","PDC API: empty results by pdc_study_id","No shared slide format (DICOM vs CZI)","Pixel sizes differ 2–5× between sources"]),
       ("What we give back",["Pre-registration with spec hashes","Pipeline code + per-slide scores","94-patient three-commons table","Public-cohort options table for GI"])]
cw,gap=3.9,0.3; x0=(13.333-(3*cw+2*gap))/2
for i,(h,items) in enumerate(heads):
    box=sl.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(x0+i*(cw+gap)),Inches(1.95),Inches(cw),Inches(4.4))
    box.fill.solid(); box.fill.fore_color.rgb=CARD; box.line.color.rgb=RGBColor(0xC3,0xD0,0xE4); box.line.width=Pt(1)
    box.shadow.inherit=False; box.text_frame.word_wrap=True
    box.text_frame.vertical_anchor=MSO_ANCHOR.TOP
    box.text_frame.margin_top=Inches(0.22); box.text_frame.margin_left=Inches(0.18); box.text_frame.margin_right=Inches(0.18)
    fill(box.text_frame,[dict(text=h,size=17,bold=True,color=NAVY,sa=12)]+
         [dict(text="•  "+it,size=14,color=DARK,sa=10) for it in items])

# ---------- 13 suggestions, deposition, next steps ----------
sl=content_slide("What we suggest, where it goes, what is next","The three items the report-out template asks for beyond results.")
heads2=[("What we suggest",["Expose MLH1 and PMS2 in the PDC quantitative matrix","Let filesPerStudy accept pdc_study_id","Publish a patient-level IDC / GDC / PDC crosswalk","Flag slides with a molecular label but no image URL"]),
        ("Where our outputs go",["Public repository: pipeline code, per-slide scores","Pre-registration with spec hashes alongside the code","Reader-study image packets to a data archive (too large for git)","All result tables released with the deck"]),
        ("What is next",["Score the remaining 301 SurGen cases: KRAS interval \u00b10.13 \u2192 \u00b10.06","Find or build an external gastric cohort \u2014 none is public today","Pathologist review of the attribution maps","MLH1 methylation as the missing sporadic-MSI leg"])]
cw,gap=3.9,0.3; x0=(13.333-(3*cw+2*gap))/2
for i,(h,items) in enumerate(heads2):
    box=sl.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(x0+i*(cw+gap)),Inches(1.95),Inches(cw),Inches(4.4))
    box.fill.solid(); box.fill.fore_color.rgb=CARD; box.line.color.rgb=RGBColor(0xC3,0xD0,0xE4); box.line.width=Pt(1)
    box.shadow.inherit=False; box.text_frame.word_wrap=True
    box.text_frame.vertical_anchor=MSO_ANCHOR.TOP
    box.text_frame.margin_top=Inches(0.22); box.text_frame.margin_left=Inches(0.18); box.text_frame.margin_right=Inches(0.18)
    fill(box.text_frame,[dict(text=h,size=16,bold=True,color=NAVY,sa=12)]+
         [dict(text="\u2022  "+it,size=13,color=DARK,sa=10) for it in items])
sl.notes_slide.notes_text_frame.text="Deposition is planned, not done: no public repository exists yet. SurGen extension is 301 more cases, about 1.2 TB of slides."

# ---------- 14 closing ----------
sl=content_slide("In one sentence")
b=sl.shapes.add_textbox(Inches(1.0),Inches(1.9),Inches(11.4),Inches(4.4))
fill(b.text_frame,[
 dict(text="The pipeline runs end to end across three NCI commons, and reproduces the 2019 benchmark.",size=24,sa=26),
 dict(text="MSI is the endpoint that travels; TP53 and KRAS do not, yet.",size=24,sa=26),
 dict(text="Code, per-slide scores and the pre-registration go out open.",size=24,sa=26),
 dict(text="Thank you — questions welcome.",size=20,color=GREY,sa=0)])
prs.save('NCI_Jamboree_ReportOut_GI_HE_molecular.pptx')
print("saved",'NCI_Jamboree_ReportOut_GI_HE_molecular.pptx',len(prs.slides.__iter__.__self__._sldIdLst))
