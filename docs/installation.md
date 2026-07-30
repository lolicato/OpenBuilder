# Installation

## 📥 Clone the repository

### Stable version

```bash
git clone https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
```

### Development version (specific branch)

```bash
git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
```

---

## 🐍 Conda Environment

We recommend installing OpenMembraneBuilder using a Conda environment:

```bash
conda env create -f environment.yml
conda activate OpenMembraneBuilder # activate OpenMembraneBuilder-dev for development version
```

---

## 🚀 Run OpenMembraneBuilder

### GUI mode

```bash
streamlit run app.py
```

### CLI mode

```bash
python app.py --no-gui outputs/OpenMembraneBuilder-xxxx/config.json
```

---

## 📦 Notes

- `coby` is used for system generation  
- `MDAnalysis` handles structure manipulation  
- `streamlit` provides the GUI  
- `streamlit-molstar` enables 3D visualization  

---

## 🔧 Optional (Documentation)

If you want to build the documentation locally:

```bash
pip install mkdocs mkdocs-material "mkdocstrings[python]"
mkdocs serve
```