from BeamLib import *

BeamStructure = Beam(1)

Shigley = BeamLab(BeamStructure)
Shigley.AddLoadCase('Shigley_11')
Shigley.LoadCases['Shigley_11'].AddFixedSupport('a', Position = 0)
Shigley.LoadCases['Shigley_11'].AddPinnedSupport('b', Position = 1)
Shigley.LoadCases['Shigley_11'].AddForce('c', Position = 0.5 ,Magnitude = -1)
Shigley.SolveAll()
Shigley.PlotAll()