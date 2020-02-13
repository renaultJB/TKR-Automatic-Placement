# TKR-Automatic-Placement

__Still under construction__

The goal of those pieces of code is to provide a way to automaticaly place the implant in the bone.
The output of the code is a placement matrix that can be used in another software to perform the placement.

For now, the method is developped for the tibia. It's based on a previous method to automatically assign anatomical coordinate systems to the tibia, femur and patella https://github.com/renaultJB/GIBOC-Knee-Coordinate-System (see [10.1016/j.jbiomech.2018.08.028](https://doi.org/10.1016/j.jbiomech.2018.08.028)). The placement is obtained with a Matlab script that account for the surgical criteria and objectives.

The tibial implant is virtually implanted on the tibia using FreeCADCmd https://www.freecadweb.org/.

*_coding languages_*
- Matlab
- Python



