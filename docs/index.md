# OpenBuilder-v1.1.0

**OpenBuilder** is a graphical and command-line tool for building Martini coarse-grained membrane and membrane–protein systems.

It supports:

- Membrane-only systems  
- Membrane systems with coarse-grained proteins  
- Martini v3 force fields  
- Reproducible JSON-based builds  
- GUI and CLI workflows  
- Automatic zipped output folders  

---

## 🚀 Main Workflow

```
GUI input → Config object → System build → config.json → ZIP archive
```

---

## 📁 Output Structure

Each build creates a self-contained project folder:

```
outputs/OpenBuilder-YYYY-MM-DD_HH-MM-SS/
├── config.json
├── user_inputs/
├── simulations/
├── toppar/
└── mdp/
```

### 📦 Zipped Output

A compressed archive is also created automatically:

```
downloads/OpenBuilder-YYYY-MM-DD_HH-MM-SS.zip
```

This archive contains everything needed to reproduce the system, including:

- Input files (`user_inputs/`)
- Configuration (`config.json`)
- Simulation setup files
