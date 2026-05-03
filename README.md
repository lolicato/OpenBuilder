# 🚀 OpenBuilder-v1.0.3-dev

![Version](https://img.shields.io/badge/version-v1.0.3--dev-red)
![Python](https://img.shields.io/badge/python-3.10-blue)
![Status](https://img.shields.io/badge/status-active-success)

🌐 **📖 Documentation:** https://lolicato.github.io/OpenBuilder  
👉 *Full user guide, tutorials, and API reference available online*

---

**OpenBuilder** is a modern **GUI + CLI toolkit for building MARTINI coarse-grained membrane and membrane–protein systems**, powered by COBY.

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

### Clone the Development version


```bash
git clone -b development https://github.com/lolicato/OpenBuilder.git
cd OpenBuilder
```

### Create environment

```bash
conda env create -f environment.yml
conda activate OpenBuilder-dev
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
python app.py --no-gui config.json
```

---

## 📁 Output Structure

Each build creates a self-contained project:

```text
outputs/OpenBuilder-xxxx/
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

👉 **Online docs (recommended):**  
https://lolicato.github.io/OpenBuilder  

👉 Local preview:

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

OpenBuilder builds upon:

- COBY  
- MDAnalysis  
- Streamlit  
- NumPy / SciPy  

