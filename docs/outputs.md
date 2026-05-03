# Output Structure

OpenBuilder writes all builds into the `outputs/` directory.

---

## 📁 Project Folder

Each build creates a self-contained project:

```
outputs/OpenBuilder-YYYY-MM-DD_HH-MM-SS/
├── config.json
├── user_inputs/
├── simulations/
├── toppar/
└── mdp/
```

---

## 📄 config.json

Stores the full configuration used for the build.

- Enables full reproducibility  
- Can be reused with CLI mode  

---

## 🧬 user_inputs/

Contains all uploaded input files:

```
user_inputs/
├── protein.pdb
└── protein.itp
```

These files are copied automatically and included in the final ZIP.

---

## 🔁 simulations/

Contains one folder per replica:

```
simulations/
├── R0001/
├── R0002/
└── R0003/
```

Each folder includes:

- generated system files  
- topology (`topol.top`)  
- simulation outputs (`.gro`, etc.)  

---

## ⚛️ toppar/

Contains force field and topology files:

- Martini force field files  
- protein topology (`protein.itp`)  

---

## 🧪 mdp/

Contains simulation parameter files:

- energy minimization (`em.mdp`)  
- equilibration (`eq*.mdp`)  
- production (`md.mdp`)  

---

## 📦 ZIP Archive

A compressed copy of the entire project is automatically created:

```
downloads/OpenBuilder-YYYY-MM-DD_HH-MM-SS.zip
```

This archive includes everything needed to reproduce the system.

---

## ✅ Summary

Each build is:

- Fully self-contained  
- Reproducible via `config.json`  
- Ready for simulation workflows  
