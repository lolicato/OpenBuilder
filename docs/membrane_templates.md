# Membrane Templates

OpenBuilder supports loading membrane compositions from **CSV** or **OpenBuilder template (`.ob`)** files. This feature allows users to quickly generate reproducible membrane systems and batch-create multiple membrane compositions from a single file.

---

## Complete Example

```text
# ExampleMembrane1

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.60,0.40,0.61,0.61
DOPC,0.25,0.40,0.67,0.67
CHOL,0.15,0.20,0.40,0.40

# ExampleMembrane2

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.60,0.20,0.61,0.61
DOPC,0.25,0.20,0.67,0.67
CHOL,0.15,0.60,0.40,0.40

```

This example defines two asymmetric membranes composed of POPC, DOPC, and cholesterol.

---



## File Structure

A template consists of one or more membrane definitions.

Each membrane definition starts with a name line:

```text
# MembraneName
```

The membrane name is used for:

- Template identification
- Replica folder naming
- Output organization

If only a single membrane composition is required, the name line can be omitted and the file may be used as a standard CSV.

---

## Template Format

Lines starting with `;` are treated as comments and ignored by OpenBuilder.
Each lipid is defined on a separate line:

```text
; resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
```

### Columns

| Column | Description |
|----------|-------------|
| `resname` | Lipid residue name |
| `upper_leaflet` | Ratio or count in the upper leaflet |
| `lower_leaflet` | Ratio or count in the lower leaflet |
| `apl_upper` | Area per lipid (nm²) in the upper leaflet |
| `apl_lower` | Area per lipid (nm²) in the lower leaflet |

---

## Single Membrane Example

```text
;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.70,0.30,0.61,0.61
DOPC,0.30,0.70,,
```

Missing APL values automatically receive the default value of **0.60 nm²**.

---

## Multiple Membrane Definitions

A single file may contain multiple membrane compositions.

```text
# PlasmaMembrane

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.70,0.30,0.61,0.61
DOPC,0.30,0.70,,

# CholesterolRich

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.50,0.40,0.61,0.61
CHOL,0.50,0.60,0.40,0.40
```

Each membrane will be processed independently.

---

## Area Per Lipid APL

APL controls membrane packing density.

| APL | Effect |
|------|--------|
| Smaller values | Tighter packing |
| Larger values | Looser packing |

If omitted:

```text
POPC,0.70,0.30
```

or

```text
POPC,0.70,0.30,,
```

OpenBuilder assigns **0.60 nm²** automatically.

---

## Ratios and Absolute Counts

Two composition modes are supported.

### Relative Ratios

```text
POPC,0.70,0.30
DOPC,0.30,0.70
```

### Absolute Counts

```text
POPC,70,30
DOPC,30,70
```

!!! warning

    Do not mix ratios and absolute counts within the same membrane definition.

---

## Supported File Types

OpenBuilder accepts:

- `.csv`
- `.ob`

Both formats use the same internal structure.

