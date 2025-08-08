# ANAV — Low-Cost 2D Position-Sensitive Detector for Muography  
**Development and Testing of a Low-Cost Two-Dimensional Detector for Muography Applications**  
**Authors:** Nitesh Kumar & Abhishek Keshri  
**Dept. of Electrical Engineering, National Institute of Technology Patna**

---

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)  
![Platform](https://img.shields.io/badge/platform-Detector%20%7C%20SiPM-green)  
![Inspiration](https://img.shields.io/badge/Inspiration-CERN%20Beamline%20for%20Schools-lightgrey)

---

## 📌 Project Summary
This repository contains code, documentation and analysis for a **compact, low-cost, two-dimensional position-sensitive detector** built from a 5×5 array of CsI scintillators read out by silicon photomultipliers (SiPMs). The detector reconstructs particle hit positions using intensity-ratio methods and is targeted at muography and particle-imaging use cases. The experimental approach and beamline validation methodology were inspired by the **CERN Beamline for Schools** program.

---

## 🔭 Table of Contents
- [Project Summary](#-project-summary)  
- [Key features](#-key-features)  
- [Repository layout](#-repository-layout)  
- [Hardware & Instrumentation](#-hardware--instrumentation)  
- [Experimental Methods](#-experimental-methods)  
- [Important Findings & Performance](#-important-findings--performance)  
- [Reproduce / Run Analysis](#-reproduce--run-analysis)  
- [Software & Dependencies](#-software--dependencies)  
- [Calibration & Analysis Recommendations](#-calibration--analysis-recommendations)  
- [Safety & Field Notes](#-safety--field-notes)  
- [Contact & License](#-contact--license)

---

## 🔬 Key features
- 5×5 CsI scintillator tile array (1 cm × 1 cm each) read by four surrounding SiPM boards.  
- Intensity-ratio position reconstruction using opposing SiPM pairs (1↔4, 2↔3).  
- Validated with long-duration cosmic-ray acquisition and controlled 2.01 GeV electron-beam scans.  
- Simple DAQ/analysis pipeline suitable for field deployment and scalable to larger arrays.

---

---

## ⚙️ Hardware & Instrumentation
- **Scintillators:** CsI tiles (1×1 cm).  
- **Photodetectors:** SiPM boards (four sides).  
- **DAQ chain:** SiPM → low-noise preamp → ADC → microcontroller/PC logger.  
- **Mechanics:** Light-tight enclosure and motorized stage for beamline scans.  
- **Test beam:** Collimated 2.01 GeV electron beam used for repeatable scans.

---

## 🧪 Experimental Methods (summary)
1. **Cosmic-ray validation** — 24-hour background run; build 2D histograms of pairwise SiPM intensity ratios to confirm positional sensitivity.  
2. **Beamline validation** — Collimated electron beam (1 cm) with motorized stage; measured position response at 1 cm steps across detector.  
3. **Reconstruction** — Ratios used:
   - `X = E2 / (E2 + E3)`  
   - `Y = E1 / (E1 + E4)`  
   Calibration maps convert ratio pairs to physical coordinates.

---

## 📈 Important Findings & Performance
- **2D imaging feasible:** Cosmic-ray histograms show distinct lines corresponding to rows/columns of the 5×5 array.  
- **One-dimensional accuracy validated:** Beam scans localized the beam to expected positions at 1 cm spacing.  
- **Non-uniform coupling:** Central tiles show stronger signals due to SiPM placement—requires per-channel calibration.  
- **Long-term stability:** Background acquisition demonstrates stable event rates after thresholding and baseline correction.  
- **Scalability:** Method scales to larger arrays with similar readout topology; calibration remains the main task.

> Full quantitative plots, calibration curves and experimental logs are in `docs/IRoC-U_2025_Qual_round_report1_Team-11318.pdf`.

---

## 🧾 How to reproduce (high level)
1. **Assemble detector:** mount 5×5 CsI array, position four SiPM boards, ensure light-tightness.  
2. **Electronics:** connect SiPM outputs to transimpedance amplifiers → ADC → microcontroller/PC.  
3. **Acquire calibration data:** dark frames, single-tile illumination or beamline hits at known positions.  
4. **Compute LUT:** map `(X_ratio, Y_ratio)` → `(x, y)`.  
5. **Run acquisition & analysis:** use scripts in `software/acquisition/` and `software/analysis/` to generate histograms and position reconstructions.

---

## 🛠 Software & Dependencies
Recommended environment:
```bash
python 3.8+
pip install numpy scipy matplotlib pandas jupyter



