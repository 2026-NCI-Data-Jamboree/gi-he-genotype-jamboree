"""Computational pathology torch stack - OpenSlide tiling, ResNet/timm models, IDC access.
torch 2.5.1 PyPI wheel bundles CUDA 12.4 -> debian_slim base."""
import modal

META = {
    "packages": ["torch", "torchvision", "timm", "openslide", "cv2", "h5py", "sklearn", "pandas", "idc_index", "transformers", "wsidicom", "imagecodecs"],
    "gpu_default": "T4",
    "egress_domains": [],
}

def build(*, secrets=None):
    img = (
        modal.Image.debian_slim(python_version="3.11")
        .pip_install("torch==2.5.1", "torchvision==0.20.1")
        .pip_install(
            "openslide-python==1.4.1", "openslide-bin",
            "timm==1.0.11", "opencv-python-headless==4.10.0.84",
            "h5py", "pandas", "scikit-learn", "pillow", "tqdm",
            "pyarrow", "idc-index", "numpy<2.2",
        )
        .pip_install("transformers==4.46.3", "huggingface_hub==0.26.5")
        .pip_install("wsidicom", "imagecodecs")
    )
    vols = {"/data": modal.Volume.from_name("idc-gi-slides", create_if_missing=True)}
    return img, vols, {"PYTHONUNBUFFERED": "1"}
