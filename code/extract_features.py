import os, io, warnings, numpy as np, pandas as pd, requests, pydicom
warnings.filterwarnings('ignore')
from pydicom.pixels import pixel_array
from PIL import Image
from skimage.color import rgb2hed
from skimage.feature import graycomatrix, graycoprops, local_binary_pattern
from skimage.morphology import remove_small_objects
from scipy.ndimage import gaussian_filter
import xml.etree.ElementTree as ET
BASE="https://idc-open-data.storage.googleapis.com/"
NS={'s':'http://doc.s3.amazonaws.com/2006-03-01'}
TARGET=4.0
def objects(uuid):
    r=requests.get(BASE,params={"prefix":uuid+"/"},timeout=120); r.raise_for_status()
    root=ET.fromstring(r.text)
    return sorted([(k.findtext('s:Key',namespaces=NS),int(k.findtext('s:Size',namespaces=NS)))
                   for k in root.findall('s:Contents',NS)], key=lambda x:x[1])
def spacing(ds):
    try: return float(ds.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing[0])*1000
    except Exception:
        try: return float(ds.PixelSpacing[0])*1000
        except Exception: return None
def header(key):
    r=requests.get(BASE+key,headers={"Range":"bytes=0-3145727"},timeout=180); r.raise_for_status()
    return pydicom.dcmread(io.BytesIO(r.content),force=True,stop_before_pixels=True)
def fetch(key,path):
    with requests.get(BASE+key,stream=True,timeout=900) as r:
        r.raise_for_status()
        with open(path,'wb') as f:
            for ch in r.iter_content(1<<20): f.write(ch)
def assemble(frames,W,H,tw,th):
    ncol=-(-W//tw); nrow=-(-H//th)
    canvas=np.full((nrow*th,ncol*tw,3),255,np.uint8)
    for i,f in enumerate(frames):
        r,c=divmod(i,ncol); canvas[r*th:(r+1)*th,c*tw:(c+1)*tw]=f
    return canvas[:H,:W]
def load_level(uuid,tmp):
    cands=[]
    for key,size in objects(uuid):
        if size < 2e5 or size > 1.2e8: continue
        try: ds=header(key)
        except Exception: continue
        sp=spacing(ds)
        if sp is None or 'THUMBNAIL' in list(getattr(ds,'ImageType',[])): continue
        if sp < 0.9: continue
        cands.append((abs(sp-TARGET), size, key, sp))
    if not cands: raise RuntimeError('no suitable level')
    cands.sort(key=lambda t:(t[0], t[1]))
    _,_,key,sp = cands[0]
    fetch(key,tmp)
    full=pydicom.dcmread(tmp)
    arr=pixel_array(full,raw=True)
    if arr.ndim==3: arr=arr[None]
    img=assemble(arr,full.TotalPixelMatrixColumns,full.TotalPixelMatrixRows,full.Columns,full.Rows)
    f=sp/TARGET
    if abs(f-1)>0.05:
        nw,nh=max(32,int(img.shape[1]*f)),max(32,int(img.shape[0]*f))
        img=np.asarray(Image.fromarray(img).resize((nw,nh),Image.BILINEAR))
    return img, sp, key
def features(img):
    f={}
    rgb=img.astype(np.float32)/255.
    gray=rgb.mean(2)
    od=-np.log10(np.clip(rgb,1/255,1).mean(2))
    mask=remove_small_objects(od>0.15,64)
    f['tissue_frac']=float(mask.mean()); f['tissue_area_mm2']=float(mask.sum()*(TARGET**2)/1e6)
    if mask.sum()<5000: raise RuntimeError('too little tissue')
    hed=rgb2hed(rgb); h,e=hed[...,0],hed[...,1]
    for nm,ch in (('hema',h),('eos',e)):
        v=ch[mask]
        f[f'{nm}_mean'],f[f'{nm}_std']=float(v.mean()),float(v.std())
        f[f'{nm}_p10'],f[f'{nm}_p50'],f[f'{nm}_p90']=[float(x) for x in np.percentile(v,[10,50,90])]
    f['hema_eos_ratio']=f['hema_mean']/(f['eos_mean']+1e-6)
    f['hema_eos_corr']=float(np.corrcoef(h[mask],e[mask])[0,1])
    f['hema_high_frac']=float((h[mask]>f['hema_mean']+f['hema_std']).mean())
    for s in (1,4): f[f'hema_localstd_s{s}']=float(gaussian_filter(h,s)[mask].std())
    for i,c in enumerate('rgb'):
        f[f'{c}_mean']=float(rgb[...,i][mask].mean()); f[f'{c}_std']=float(rgb[...,i][mask].std())
    mx,mn=rgb.max(2),rgb.min(2); sat=np.where(mx>0,(mx-mn)/np.maximum(mx,1e-6),0)
    f['sat_mean'],f['sat_std']=float(sat[mask].mean()),float(sat[mask].std())
    f['val_mean'],f['val_std']=float(mx[mask].mean()),float(mx[mask].std())
    g8=np.where(mask,np.digitize(gray,np.linspace(0,1,9)[1:-1]).astype(np.uint8),0)
    gl=graycomatrix(g8,[1,2,4,8],[0,np.pi/4,np.pi/2,3*np.pi/4],levels=9,symmetric=True,normed=True)
    gl=gl[1:,1:]; gl=gl/np.maximum(gl.sum((0,1),keepdims=True),1e-9)
    for prop in ('contrast','homogeneity','correlation','ASM','dissimilarity'):
        vals=graycoprops(gl,prop)
        for di,dd in enumerate([1,2,4,8]): f[f'glcm_{prop}_d{dd}']=float(vals[di].mean())
    lbp=local_binary_pattern((gray*255).astype(np.uint8),8,1,'uniform')
    hist=np.bincount(lbp[mask].astype(int),minlength=10)[:10].astype(float); hist/=hist.sum()
    for i,v in enumerate(hist): f[f'lbp_{i}']=float(v)
    return f
def one(rec):
    uuid=rec['crdc_series_uuid']; tmp=f"/tmp/{uuid}.dcm"
    try:
        img,sp,key=load_level(uuid,tmp)
        f=features(img); f.update(rec); f['level_um_per_px']=sp; f['ok']=True; f['error']=''
    except Exception as ex:
        f=dict(rec); f['ok']=False; f['error']=f"{type(ex).__name__}: {ex}"
    finally:
        if os.path.exists(tmp): os.remove(tmp)
    return f
if __name__=='__main__':
    recs=pd.read_csv('pilot_slides.csv').to_dict('records')
    done=set()
    if os.path.exists('pilot_features.csv'):
        prev=pd.read_csv('pilot_features.csv'); prev=prev[prev.ok==True]
        done=set(prev.crdc_series_uuid); out=prev.to_dict('records')
    else: out=[]
    recs=[r for r in recs if r['crdc_series_uuid'] not in done]
    print('to do:',len(recs),flush=True)
    from multiprocessing import Pool
    with Pool(4) as p:
        for i,r in enumerate(p.imap_unordered(one,recs),1):
            out.append(r)
            if i%10==0:
                pd.DataFrame(out).to_csv('pilot_features.csv',index=False)
                print(i,'of',len(recs),'| failed',sum(1 for x in out if not x['ok']),flush=True)
    pd.DataFrame(out).to_csv('pilot_features.csv',index=False)
    print('FINISHED rows',len(out),'failed',sum(1 for x in out if not x['ok']),flush=True)
