# 🚀 OpenMembraneBuilder

![Version](https://img.shields.io/badge/version-v1.0.4--dev-red)
![Python](https://img.shields.io/badge/python-3.10-blue)
![Status](https://img.shields.io/badge/status-active-success)

**OpenMembraneBuilder** is a modern **GUI + CLI toolkit for building MARTINI coarse-grained membrane and membrane–protein systems**, powered by COBY.

🌐 **Documentation:** https://lolicato.github.io/OpenMembraneBuilder

It combines:
- 🧬 **MDAnalysis** for structure handling  
- ⚙️ **COBY** for system generation  
- 🌐 **Streamlit** for an interactive GUI  

with a strong focus on **reproducibility**, **usability**, and **local execution**.

---

## 🖼️ Interface

![Application GUI](./pictures/GUI.png)

---

## ✨ Features

- 🧪 Build membrane and membrane–protein systems  
- 🎛️ Fully interactive Streamlit GUI  
- 🧾 Automatic **JSON-based reproducibility**  
- 💻 Command-line interface (CLI) for batch workflows  
- 📦 Self-contained outputs with input files included  
- 🔁 Multi-replica system generation  

---

## 🧭 Workflow

```text
GUI → Config → Build → config.json → ZIP → CLI reuse
```

---

## 📦 Installation

### Clone the repository

```bash
git clone https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
```

### (Optional) Development branch

```bash
git clone -b dev https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
```

### Create environment

```bash
conda env create -f environment.yml
conda activate OpenMembraneBuilder-dev
```

---

## 🚀 Usage

### GUI mode

```bash
streamlit run app.py
```

---

### CLI mode

```bash
python app.py --no-gui outputs/OpenMembraneBuilder-xxxx/config.json
```

---

## 📁 Output Structure

Each build creates a self-contained project:

```text
outputs/OpenMembraneBuilder-xxxx/
├── config.json
├── user_inputs/
├── simulations/
├── toppar/
└── mdp/
```

✔ Fully reproducible  
✔ Ready for simulation  
✔ ZIP archive generated automatically  

---

## 🔁 Reproducibility

Every build generates:

```text
config.json
```

This allows:

- exact rebuilds  
- parameter editing  
- automated workflows  

---

## 📚 Documentation

Full documentation available online:

👉 https://lolicato.github.io/OpenMembraneBuilder

Or run locally:

```bash
mkdocs serve
```

---

## 🤝 Contributing

We welcome:

- Bug reports  
- Feature requests  
- Documentation improvements  
- New workflows  

---

## 📄 License

Currently restricted to **academic and personal use**.

Open-source release planned after first publication.

---

## 🙏 Acknowledgements

OpenMembraneBuilder builds upon:

- COBY  
- MDAnalysis  
- Streamlit  
- NumPy / SciPy  

---

## 🌐 Future

- Hosted version: **open-builder.com**  
- Advanced analysis workflows  
- Multi-system automation  
- Improved visualization  

---

## ⭐ If you like this project

Consider starring the repository!
