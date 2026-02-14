# H2: Adaptive Kernel Framework for Galaxy Rotation Curves

**Status:** In development (Phase-4 data integrity issue being resolved)

**Background:** This project tests whether adaptive nonlocal kernel 
deformation can improve inner-region rotation curve diagnostics 
identified in previous work (H1.5).

**Methodology:** Pre-registered diagnostic tests (see Diaries/)

**Current Issue:** Basis grid regeneration required due to 
mass-to-light ratio mismatch (see diagnostics/phase4_checkpoint.md)

**Contact:** [your email]
```

---

## 🎯 **GitHub Strategy:**

### **Repository Name:**
`H2-adaptive-kernels` (descriptive, searchable, no Ananta reference)

### **License:**
MIT or GPL-3 (allows others to use/cite but requires attribution)

### **Repository Structure:**
```
H2-adaptive-kernels/
├── README.md (project status, clear "in development" label)
├── diaries/ (all 7 H2 Logic Diaries)
├── core/ (galaxy_io, h2_engine, kernels, etc)
├── diagnostics/ (all test scripts)
├── data/
│   ├── sparc/ (SPARC input data)
│   └── derived/ (outputs, including "broken" Phase-4)
├── docs/
│   └── phase4_checkpoint.md (document current data issue)
└── LICENSE
```

### **First Commit Message:**
```
Initial commit: H2 framework with pre-registered diagnostics

Phase-4 in progress. Current blocker: basis grid parameter mismatch
identified in Test-2 correlation diagnostic. See docs/phase4_checkpoint.md
