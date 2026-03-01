# set_wireframe.py
from paraview.simple import *
view = GetActiveViewOrCreate('RenderView')
# ensure all currently shown sources use Wireframe
for rep in GetRepresentations():
    rep.Representation = 'Wireframe'
Render()

