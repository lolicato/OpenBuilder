# Installation

## 📥 Clone the repository

### Stable version

```bash
git clone https://github.com/lolicato/OpenBuilder.git
cd OpenBuilder
```

### Development version (specific branch)

```bash
git clone -b development https://github.com/lolicato/OpenBuilder.git
cd OpenBuilder
```

---

## 🐍 Conda Environment

We recommend installing OpenBuilder using a Conda environment:

```bash
conda env create -f environment.yml
conda activate OpenBuilder # activate OpenBuilder-dev for development version
```

---

## 🚀 Run OpenBuilder

### GUI mode

```bash
streamlit run app.py
```

### CLI mode

```bash
python app.py --no-gui outputs/OpenBuilder-xxxx/config.json
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
