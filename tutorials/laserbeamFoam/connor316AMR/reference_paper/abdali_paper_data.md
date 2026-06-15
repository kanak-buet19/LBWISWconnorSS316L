import os

# Create the updated content for the markdown file including geometry and mesh size
md_content_v2 = """# Selective Laser Melting (SLM) Solidification Data

## 1. Study Reference & Authors
* **Authors:** Amirreza Abdali, Syamak Hossein Nedjad, Habib Hamed Zargari, Abdollah Saboori, Mehmet Yildiz
* **Journal:** Journal of Materials Research and Technology (2024)

## 2. Material Specifications
* **Base Material:** Hot-rolled AISI 316L Stainless Steel plate
* **Plate Thickness:** 12 mm
* **Chemical Composition (wt.%):** 0.019C - 0.514Si - 1.32Mn - 16.51Cr - 10.46Ni - 2.14Mo - 0.401Cu - 0.019Ti - 0.089V - 0.025P - 0.008S - 68.50Fe

## 3. Simulation Geometry & Mesh Size (SLM Model)
* **Computational Domain Size (Geometry):** 1.2 mm × 0.3 mm × 0.3 mm
* **Mesh Size (Grid Resolution):** Discretized using rectangular cuboid cells of **5 $\\mu$m**

## 4. Process Parameter Set (SLM Case for Figure 5b)
The localized cooling rate depth profile reported in the study corresponds to the following specific operational parameter set:
* **Laser Power ($P$):** 95 W
* **Laser Beam Spot Size (Radius):** 40 $\\mu$m
* **Scan Speed ($V$):** 200 mm/s
* **Substrate Configuration:** Powderless / Bare plate (Autogenous laser track)

## 5. Individual Localized Cooling Rate Values (Figure 5b)
The individual non-averaged cooling rates computed at specific distances from the melt pool boundary (bottom) are listed below:

| Distance from Melt Pool Bottom ($\\mu$m) | Computed Solidification Cooling Rate (K/s) |
| :---: | :---: |
| $\\sim$ 20 | 972,400 |
| $\\sim$ 80 | 384,300 |
| $\\sim$ 140 | 319,500 |
| $\\mu$m $\\sim$ 200 | 286,400 |

*Note: These local values represent non-averaged points extracted via 3D finite difference numerical simulations mapping the thermal profile between the liquidus (1700 K) and solidus (1650 K) isotherms.*
"""

filename_v2 = "slm_cooling_rate_data-v2.md"
with open(filename_v2, "w", encoding="utf-8") as f:
    f.write(md_content_v2)

print(f"File generated successfully: {filename_v2}")