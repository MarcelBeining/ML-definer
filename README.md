# ML-definer

This small graphical user interface helps to map the anatomical borders of dentate gyrus in an imaged volume onto morphologies that had been reconstructed in that volume.
It is one of the tools used (and should be cited like this) for the paper [Beining M*, Jungenitz T*, Radic T, Deller T, Cuntz H, Jedlicka P, Schwarzacher SW: Adult-born dentate granule cells show a critical period of dendritic reorganization and are distinct from developmentally born cells. Brain Structure and Function 2017, 222(3):1427-1446.](https://link.springer.com/article/10.1007/s00429-016-1285-y)

The tool allows to draw the two borders of the granule cell layer as well as the hippocampal fissure in several depth of the imaged volume. Borders between the outer, middle, and inner ML (OML, MML, and IML, respectively) 
are then automatically calculated by dividing the area between the fissure and the GCL into three equal parts. Also the borders are inter/extrapolated along the z dimension.
Finally nodes of reconstructed trees are automatically assigned to one of the following regions according to the layer contours: subgranular zone (SGZ), GCL, IML, MML, or OML.

# Getting started
Start the gui_getMLcontours_big.m file in Matlab and follow the instructions.

**Keys:**
- Left or right mouse key button lets you draw a point
- F and D lets you move forward and backward along z in the image volume.
- Backspace lets you remove the last point
- Space lets you finish the current drawing

# Licensing
This software is published under the MIT license. For further information, read the LICENSE file.
