# H2: Adaptive Kernel Framework for Galaxy Rotation Curves

**Status:** In development (Phase-4 data integrity issue being resolved)

**Background:** This project tests whether adaptive nonlocal kernel 
deformation can improve inner-region rotation curve diagnostics 
identified in previous work (H1.5).

**Methodology:** Pre-registered diagnostic tests (see Diaries/)

**Current Issue:** Basis grid regeneration required due to 
mass-to-light ratio mismatch (see diagnostics/phase4_checkpoint.md)

**Contact:** [vitor.figueiredo.research@protonmail.com]
```

---

H2-adaptive-kernels

**License:**
MIT

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


