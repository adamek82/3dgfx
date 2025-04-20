# 3D Graphics Experiments

*Simple experiments with 3D graphics using plain JavaScript and HTML5 Canvas*

---

## Overview

Two interactive browser demos showcasing different 3D→2D projection techniques:

- **blocks.html**  
  A “glass skyscraper” point‑cloud demo using a basic projection algorithm.

- **apka.html**  
  A full 3D engine with `Point3d`/`Vector`/`Matrix` math, observer coordinate systems, virtual screen projection, and a rendered cube.

No build step required—just open in your browser and play with the sliders to change camera and canvas parameters in real time.

---

## Features

- **Basic projection** (`3dfx.js`)  
  - Computes eye‑space coordinates and projects 3D points onto 2D canvas  
  - Randomized skyscraper walls demo

- **Advanced 3D engine** (`3dfx_4_app.js`)  
  - `Point2d` / `Point3d` / `Vector` / `Matrix` / `Scale` classes  
  - `C3DView` for full 3D viewing transformations, inverse mapping, and virtual screens  
  - Renders coordinate axes, observer frustum, observer’s screen, and a 3D cube

- **Interactive HTML**  
  - `blocks.html`: sliders for canvas size and observer (X, Y, Z)  
  - `apka.html`: range inputs for R, φ, θ, screen distance, observer parameters

- **Styling** (`mstyle.css`)  
  - Responsive layout, canvas container decoration, and control panel styling

---

## Project Structure

```
/3dgfx
│
├── blocks.html          # Basic demo: random point‑cloud skyscraper
├── apka.html            # Advanced demo: cube, axes, virtual screen
├── mstyle.css           # Shared CSS for layout & styling
│
├── 3dfx.js              # Core projection & skyscraper logic
├── 3dfx_4_app.js        # Full 3D math library & rendering engine
│
└── 3dgfx.code-workspace # VS Code workspace settings (optional)
```

---

## Getting Started

1. **Clone or download** this repository.  
2. **Open** `blocks.html` or `apka.html` in any modern browser (Chrome, Firefox, Edge).  
3. **Adjust** the on‑page sliders to tweak canvas dimensions and camera/observer parameters.

---

## Dependencies

- [jQuery 3.x](https://jquery.com/) (loaded via CDN in the HTML)

---

## License

This project is provided “as‑is” for learning and experimentation. Feel free to fork and extend!
