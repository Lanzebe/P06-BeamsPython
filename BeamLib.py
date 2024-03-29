from math import radians
import numpy as np
import matplotlib.pyplot as plt

class BeamLab():
    def __init__(self,Beam) -> None:
        self.Beam = Beam
        self.LoadCases = {}

        ### Simulation Settings:
        self.MinimumNumElementsForPlot = 1000

    def AddLoadCase(self,LoadCaseName):
        self.LoadCases.update({LoadCaseName: LoadCase(LoadCaseName,self)})

    def SolveAll(self):
        for LC in self.LoadCases:
            self.LoadCases[LC].SolveLoadCase()

    def PlotAll(self):
        self.xs = np.linspace(self.Beam.X_Start,self.Beam.X_End,self.MinimumNumElementsForPlot)

        distlod = []
        shearforce = []
        bendingmoment = []
        theta_angle = []
        displacement = []
        

        for LC in self.LoadCases: 
            distlod.append(self.LoadCases[LC].CalcDistributedLoading())
            shearforce.append(self.LoadCases[LC].CalcShearForce())
            bendingmoment.append(self.LoadCases[LC].CalcBendingMoment())
            theta_angle.append(self.LoadCases[LC].CalcAngleTheta())
            displacement.append(self.LoadCases[LC].CalcDisplacement())

        labels = list(self.LoadCases.keys())

        ### Dist Loading
        plt.figure("Distributed Loading")
        plt.xlabel('Beam Lenght [m]')
        plt.ylabel('Distributed Loading [N/m]')
        plt.plot(self.xs, np.transpose(distlod), label = labels)
        plt.plot([self.Beam.X_Start,self.Beam.X_End],[0,0],'k')
        plt.grid()
        plt.xlim([self.Beam.X_Start,self.Beam.X_End])
        plt.legend()
        #plt.show()

        ### Shear Force Diagram
        plt.figure("Shear Force Diagram")
        plt.xlabel('Beam Lenght [m]')
        plt.ylabel('Shear Force [N]')
        plt.plot(self.xs,np.transpose(shearforce), label = labels)
        plt.plot([self.Beam.X_Start,self.Beam.X_End],[0,0],'k')
        plt.grid()
        plt.xlim([self.Beam.X_Start,self.Beam.X_End])
        plt.legend()
        #plt.show()

        ### Bending Moments Diagram
        plt.figure("Bending Moments Diagram")
        plt.xlabel('Beam Lenght [m]')
        plt.ylabel('Bending Moment [Nm]')
        plt.plot(self.xs,np.transpose(bendingmoment), label = labels)
        plt.plot([self.Beam.X_Start,self.Beam.X_End],[0,0],'k')
        plt.grid()
        plt.xlim([self.Beam.X_Start,self.Beam.X_End])
        plt.legend()
        #plt.show()

        ### Theta Diagram
        plt.figure("Angle of Beam")
        plt.xlabel('Beam Lenght [m]')
        plt.ylabel('Theta angle [deg]')
        plt.plot(self.xs,np.transpose(theta_angle), label = labels)
        plt.plot([self.Beam.X_Start,self.Beam.X_End],[0,0],'k')
        plt.grid()
        plt.xlim([self.Beam.X_Start,self.Beam.X_End])
        plt.legend()
        #plt.show()

        ### Displacment Diagram
        plt.figure("Displacement of Beam")
        plt.xlabel('Beam Lenght [m]')
        plt.ylabel('Deflection [m]')
        plt.plot(self.xs,np.transpose(displacement), label = labels)
        plt.plot([self.Beam.X_Start,self.Beam.X_End],[0,0],'k')
        plt.grid()
        plt.xlim([self.Beam.X_Start,self.Beam.X_End])
        plt.legend()
        plt.show()

class Beam():
    def __init__(self,X_End ,X_Start = 0,ElasticModulus = 210E9,I_Value= (4E6)*1E-12) -> None:  
        self.Lenght = X_End - X_Start
        self.X_End = X_End
        self.X_Start = X_Start
        self.ElasticModulus = ElasticModulus    #[Pa]
        self.I_Value = [[X_Start,X_End],[I_Value,I_Value]] #[m^4]

        #print('ElasticModulus:', ElasticModulus)


    def ReturnIvalue(self,X_pos):
        ### Returns the active Element for a given point
        return np.interp(X_pos,self.I_Value[0],self.I_Value[1])

class LoadCase():
    def __init__(self,LoadCaseName,LabEnvironment) -> None:
        self.LoadCaseName = LoadCaseName
        self.KeyPositions = []

        self.DistributedLoadings = {}
        self.Forces = {}
        self.Moments = {}
        self.FixedSupports = {}
        self.PinnedSupports = {}
        self.VerticalSliderSupports = {}

        

        self.BoundaryConditions = [SumOfForcesBoundaryCondition(self),SumOfMomentsBoundaryCondition(self)]

        self.Unknowns = []
        self.Lab = LabEnvironment


        self.C1 = UnknownConstant('C1',self)
        self.C2 = UnknownConstant('C2',self)
        
    ##### Building of Loadcase:
    
    def AddDistributedLoading(self,LoadingName,Magnitudes,Positions):
        self.DistributedLoadings.update({LoadingName: DistributedLoading(LoadingName,Magnitudes,Positions,self)})

    def AddForce(self,ForceName,Magnitude,Position):
        self.Forces.update({ForceName: Force(ForceName,Magnitude,Position,self)})

    def AddMoment(self,MomentName,Magnitude,Position):
        self.Moments.update({MomentName: Moment(MomentName,Magnitude,Position,self)})

    def AddFixedSupport(self,FixedSupportName,Position):
        self.FixedSupports.update({FixedSupportName: FixedSupport(FixedSupportName,Position,self)})
    
    def AddPinnedSupport(self,PinnedSupportName,Position):
        self.PinnedSupports.update({PinnedSupportName: PinnedSupport(PinnedSupportName,Position,self)})
    
    def AddVerticalSliderSupport(self,VerticalSliderSupportName,Position):
        self.VerticalSliderSupports.update({VerticalSliderSupportName: VerticalSliderSupport(VerticalSliderSupportName,Position,self)})

        ### Discretization

    def InitializeElements(self):
        for pos in self.Lab.Beam.I_Value[0]:
            self.KeyPositions.append(pos)

        self.KeyPositions = [*set(self.KeyPositions)]#.sort()
        self.KeyPositions.sort()

        self.Elements = []

        for idx, keypos in enumerate(self.KeyPositions):
            if idx !=0:
                self.Elements.append(BeamElement(self.KeyPositions[idx-1],self.KeyPositions[idx], self))

        ### Linking the Elements with boundary conditions
        for idx in range(len(self.Elements)):
            if idx !=0:
                #self.BoundaryConditions.append(LinkAngleElementBoundaryCondition(self.Elements[idx-1],self.Elements[idx],self))
                #self.BoundaryConditions.append(LinkDisplacementElementBoundaryCondition(self.Elements[idx-1],self.Elements[idx],self))
                pass

        N_Dimension = len(self.Unknowns)
        self.Matrix = np.zeros((N_Dimension,N_Dimension))
        self.Vector = np.zeros(N_Dimension)

        ### Initialize Boundary conditions:
    def InitializeBoundaryConditions(self):

        for BCdx in range(len(self.BoundaryConditions)):
            #print(self.BoundaryConditions[BCdx])
            self.Matrix[BCdx], self.Vector[BCdx] = self.BoundaryConditions[BCdx].EvaluateEquation()

    def CalcDistributedLoading(self):
        dist_lod = np.zeros(self.Lab.MinimumNumElementsForPlot)

        for element in self.Elements:
            new_xs = self.Lab.xs.copy()
            new_xs[new_xs<element.StartPos] = None
            new_xs[new_xs>=element.EndPos] = None

            dist_lod = dist_lod + np.nan_to_num(element.DistributedLoadingFunction(new_xs),0)

        return dist_lod

    def CalcShearForce(self):
        ShearForce = np.zeros(self.Lab.MinimumNumElementsForPlot)

        for element in self.Elements:
            new_xs = self.Lab.xs.copy()
            new_xs[new_xs<element.StartPos] = None
            new_xs[new_xs>=element.EndPos] = None

            ShearForce = ShearForce + np.nan_to_num(element.ShearForceFunction(new_xs),0)

        return ShearForce

    def CalcBendingMoment(self):
        BendingMoment = np.zeros(self.Lab.MinimumNumElementsForPlot)

        for element in self.Elements:
            new_xs = self.Lab.xs.copy()
            new_xs[new_xs<element.StartPos] = None
            new_xs[new_xs>=element.EndPos] = None

            BendingMoment = BendingMoment + np.nan_to_num(element.BendingMomentFunction(new_xs),0)

        return BendingMoment

    def CalcAngleTheta(self):
        AngleTheta = np.zeros(self.Lab.MinimumNumElementsForPlot)

        for element in self.Elements:
            new_xs = self.Lab.xs.copy()
            new_xs[new_xs<element.StartPos] = None
            new_xs[new_xs>=element.EndPos] = None

            AngleTheta = AngleTheta + np.nan_to_num(element.AngleThetaFunction(new_xs),0)

        return AngleTheta

    def CalcDisplacement(self):
        Displacement = np.zeros(self.Lab.MinimumNumElementsForPlot)

        for element in self.Elements:
            new_xs = self.Lab.xs.copy()
            new_xs[new_xs<element.StartPos] = None
            new_xs[new_xs>=element.EndPos] = None

            Displacement = Displacement + np.nan_to_num(element.DisplacementFunction(new_xs),0)

        return Displacement

    def SolveLoadCase(self):
        self.InitializeElements()
        self.InitializeBoundaryConditions()
        
        self.SolVector = np.linalg.solve(self.Matrix,self.Vector)

        print(self.LoadCaseName," :")
        for idx, unknown in enumerate(self.Unknowns):
            unknown.Magnitude = self.SolVector[idx]

        for idx, unknown in enumerate(self.Unknowns[2:]):
            print(unknown.Name,'  :', unknown.Magnitude)

    def ReturnElement(self,X_pos):
        ### Returns the active Element for a given point
        
        for Element in self.Elements:
            if Element != self.Elements[-1]:
                if (Element.StartPos<= X_pos) and (Element.EndPos > X_pos):
                    return Element
            else:
                if (Element.StartPos<= X_pos) and (Element.EndPos >= X_pos):
                    return Element

    def ReturnIValue(self,X_pos):
        ### Returns the Ivalue at a given points
        for idx in range(len(self.Lab.Beam.I_Value[0])):
            if idx != 0:
                if (self.Lab.Beam.I_Value[0][idx-1] <= X_pos) and (self.Lab.Beam.I_Value[0][idx] >= X_pos):
                    return self.Lab.Beam.I_Value[1][idx-1]+((X_pos-self.Lab.Beam.I_Value[0][idx-1])/(self.Lab.Beam.I_Value[0][idx]-self.Lab.Beam.I_Value[0][idx-1]))*(self.Lab.Beam.I_Value[1][idx]-self.Lab.Beam.I_Value[1][idx-1])
             
    def __str__(self) -> str:
        print('###### '+ self.LoadCaseName + ' ######')
        print('Has '+ str(len(self.Forces))+ ' force/s:')
        for force in self.Forces:
            print('  -',force)
        print('Has '+ str(len(self.Moments))+ ' Moment/s:')
        for moment in self.Moments:
            print('  -',moment)
        print('Has '+ str(len(self.FixedSupports))+ ' Fixed Support/s:')
        for fixed_support in self.FixedSupports:
            print('  -',fixed_support)
        print('Has '+ str(len(self.PinnedSupports))+ ' Pinned Support/s:')
        for pinned_support in self.PinnedSupports:
            print('  -',pinned_support)
        print('Has '+ str(len(self.VerticalSliderSupports))+ ' Vertical Slider Support/s:')
        for vertical_slider_support in self.VerticalSliderSupports:
            print('  -',vertical_slider_support)
        return'#################'

class DistributedLoading():
    def __init__(self,LoadingName,Magnitudes,Positions,LoadCaseEnvironment) -> None:
        self.Name = LoadingName
        self.Magnitudes = Magnitudes
        self.Positions = Positions
        self.LoadCase = LoadCaseEnvironment

        

        ### Adding Element Location
        for i in self.Positions:
            self.LoadCase.KeyPositions.append(i)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitudes) + '[Nm] at '+ str(self.Positions) + '[m]'

class Force():
    def __init__(self,ForceName,Magnitude,Position,LoadCaseEnvironment) -> None:
        self.Name = ForceName
        self.Magnitude = Magnitude
        self.Position = Position
        self.LoadCase = LoadCaseEnvironment
        self.CheckIfOnBeam()
        ### Adding Element Location
        self.LoadCase.Forces.update({ForceName: self})
        self.LoadCase.KeyPositions.append(Position)

    def CheckIfOnBeam(self):
        if not ((self.Position>= self.LoadCase.Lab.Beam.X_Start) and (self.Position<= self.LoadCase.Lab.Beam.X_End)):
            print('Invalid Position Argument. Place ' + self.Name + ' on the Beam, Between', self.LoadCase.Lab.Beam.X_Start, 'and', self.LoadCase.Lab.Beam.X_End)


    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[Nm] at '+ str(self.Position) + '[m]'

class Moment():
    def __init__(self,MomentName,Magnitude,Position,LoadCaseEnvironment) -> None:
        self.Name = MomentName
        self.Magnitude = Magnitude
        self.Position = Position
        self.LoadCase = LoadCaseEnvironment
        self.CheckIfOnBeam()

        ### Adding Element Location
        self.LoadCase.KeyPositions.append(Position)
        self.LoadCase.Moments.update({MomentName: self})

    def CheckIfOnBeam(self):
        if not ((self.Position>= self.LoadCase.Lab.Beam.X_Start) and (self.Position<= self.LoadCase.Lab.Beam.X_End)):
            print('Invalid Position Argument. Place ' + self.Name + ' on the Beam, Between', self.LoadCase.Lab.Beam.X_Start, 'and', self.LoadCase.Lab.Beam.X_End)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[N] at '+ str(self.Position) + '[m]'

class ReactionForce():
    def __init__(self,ForceName,Position,LoadCase) -> None:
        self.Magnitude = UnknownConstant(ForceName,LoadCase)
        self.Name = ForceName
        self.Position = Position
        self.LoadCase = LoadCase
        self.CheckIfOnBeam()
        ### Adding Element Location
        self.LoadCase.KeyPositions.append(Position)
        self.LoadCase.Forces.update({ForceName: self})
        #self.LoadCase.Unknowns.append(self)

        #self.VariableNumber = len(self.LoadCase.Unknowns)

    def CheckIfOnBeam(self):
        if not ((self.Position>= self.LoadCase.Lab.Beam.X_Start) and (self.Position<= self.LoadCase.Lab.Beam.X_End)):
            print('Invalid Position Argument. Place ' + self.Name + ' on the Beam, Between', self.LoadCase.Lab.Beam.X_Start, 'and', self.LoadCase.Lab.Beam.X_End)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[Nm] at '+ str(self.Position) + '[m]'

class ReactionMoment():
    def __init__(self,MomentName,Position,LoadCase) -> None:
        self.Magnitude = UnknownConstant(MomentName,LoadCase)
        self.Name = MomentName
        self.Position = Position
        self.LoadCase = LoadCase
        self.CheckIfOnBeam()

        ### Adding Element Location
        self.LoadCase.KeyPositions.append(Position)
        #self.LoadCase.Unknowns.append(self)
        self.LoadCase.Moments.update({MomentName: self})

        #self.VariableNumber = len(self.LoadCase.Unknowns)

    def CheckIfOnBeam(self):
        if not ((self.Position>= self.LoadCase.Lab.Beam.X_Start) and (self.Position<= self.LoadCase.Lab.Beam.X_End)):
            print('Invalid Position Argument. Place ' + self.Name + ' on the Beam, Between', self.LoadCase.Lab.Beam.X_Start, 'and', self.LoadCase.Lab.Beam.X_End)


    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[N] at '+ str(self.Position) + '[m]'

class FixedSupport():
    def __init__(self,FixedSupportName,Position,Loadcase) -> None:
        self.Name = FixedSupportName
        self.Loadcase = Loadcase
        self.Force = ReactionForce('F_' + FixedSupportName ,Position, self.Loadcase)
        self.Moment = ReactionMoment('M_' + FixedSupportName ,Position, self.Loadcase)
        self.Position = Position

        self.Loadcase.BoundaryConditions.append(DisplacementBoundaryCondition(Position,self.Loadcase))
        self.Loadcase.BoundaryConditions.append(AngleBoundaryCondition(Position,self.Loadcase))

    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Force)
        print('    -',self.Moment)
        return '_______________'

class PinnedSupport():
    def __init__(self,PinnedSupportName,Position,Loadcase) -> None:
        self.Name = PinnedSupportName
        self.Loadcase = Loadcase
        self.Force = ReactionForce('F_' +PinnedSupportName ,Position,self.Loadcase)
        self.Position = Position

        self.Loadcase.BoundaryConditions.append(DisplacementBoundaryCondition(Position,self.Loadcase))
        


    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Force)
        return '_______________'

class VerticalSliderSupport():
    def __init__(self,VerticalSliderSupportName,Position,Loadcase) -> None:
        self.Name = VerticalSliderSupportName
        self.Loadcase = Loadcase
        self.Moment = ReactionMoment('M_' +VerticalSliderSupportName,Position,self.Loadcase)
        self.Position = Position

        self.Loadcase.BoundaryConditions.append(AngleBoundaryCondition(Position,self.Loadcase))


    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Moment)
        return '_______________'

class ShearForceBoundaryCondition():
    def __init__(self,Position, LoadCase, ShearForce = 0) -> None:
        self.Position = Position
        self.Value = ShearForce
        self.LoadCase = LoadCase
        self.LoadCase.KeyPositions.append(Position)

class BendingMomentBoundaryCondition():
    def __init__(self,Position, LoadCase, BendingMoment = 0) -> None:
        self.Position = Position
        self.Value = BendingMoment
        self.LoadCase = LoadCase
        self.LoadCase.KeyPositions.append(Position)
        
class AngleBoundaryCondition():
    def __init__(self,Position,LoadCase, AngleDeg = 0) -> None:
        self.Position = Position
        self.Value = radians(AngleDeg)
        self.LoadCase = LoadCase
        self.LoadCase.KeyPositions.append(Position)

    def EvaluateEquation(self):
        ActiveElement = self.LoadCase.ReturnElement(self.Position)
        return ActiveElement.EvaluateAngleEquation(self.Position)

class DisplacementBoundaryCondition():
    def __init__(self,Position, LoadCase, Displacement = 0) -> None:
        self.Position = Position
        self.Value = Displacement
        self.LoadCase = LoadCase
        self.LoadCase.KeyPositions.append(Position)

    def EvaluateEquation(self):
        ActiveElement = self.LoadCase.ReturnElement(self.Position)

        return ActiveElement.EvaluateDeflectionEquation(self.Position)

class LinkAngleElementBoundaryCondition():
    def __init__(self, Element1, Element2, LoadCase) -> None:
        self.Element1 = Element1
        self.Element2 = Element2
        self.LoadCase = LoadCase

    def EvaluateEquation(self):
        LHS1,RHS1 = self.Element1.EvaluateAngleEquation(1)
        LHS2,RHS2 = self.Element2.EvaluateAngleEquation(0)

        return LHS2-LHS1, RHS2-RHS1

class LinkDisplacementElementBoundaryCondition():
    def __init__(self, Element1, Element2, LoadCase) -> None:
        self.Element1 = Element1
        self.Element2 = Element2
        self.LoadCase = LoadCase

    def EvaluateEquation(self):
        LHS1,RHS1 = self.Element1.EvaluateDeflectionEquation(1)
        LHS2,RHS2 = self.Element2.EvaluateDeflectionEquation(0)

        return LHS2-LHS1, RHS2-RHS1

class SumOfForcesBoundaryCondition():
    def __init__(self, LoadCase) -> None:
        self.LoadCase = LoadCase

    def EvaluateEquation(self):
        RHS = 0
        LHS = np.zeros(len(self.LoadCase.Unknowns))

        for i in self.LoadCase.Elements:
            RHS = RHS - (i.CumulativeShearForce)

        #for dist in self.LoadCase.DistributedLoadings:
        #    for posdx in range(len(self.LoadCase.DistributedLoadings[dist].Positions)):
        #        if posdx!=0:
        #            RHS = RHS - np.mean([self.LoadCase.DistributedLoadings[dist].Magnitudes[posdx-1],self.LoadCase.DistributedLoadings[dist].Magnitudes[posdx]])*(self.LoadCase.DistributedLoadings[dist].Positions[posdx]-self.LoadCase.DistributedLoadings[dist].Positions[posdx-1])

        for force in self.LoadCase.Forces:
            if type(self.LoadCase.Forces[force]) == Force:
                RHS = RHS - self.LoadCase.Forces[force].Magnitude
            
            if type(self.LoadCase.Forces[force]) == ReactionForce:
                LHS[self.LoadCase.Forces[force].Magnitude.VariableNumber-1] = 1

        return LHS, RHS

class SumOfMomentsBoundaryCondition():
    def __init__(self, LoadCase) -> None:
        self.LoadCase = LoadCase

    def EvaluateEquation(self):
        RHS = 0
        LHS = np.zeros(len(self.LoadCase.Unknowns))

        for i in self.LoadCase.Elements:
            RHS = RHS - (i.CGx*i.CumulativeShearForce)

        for force in self.LoadCase.Forces:
            if type(self.LoadCase.Forces[force]) == Force:
                RHS = RHS - self.LoadCase.Forces[force].Magnitude*self.LoadCase.Forces[force].Position
            
            if type(self.LoadCase.Forces[force]) == ReactionForce:
                LHS[self.LoadCase.Forces[force].Magnitude.VariableNumber-1] = self.LoadCase.Forces[force].Position

        for moment in self.LoadCase.Moments:
            if type(self.LoadCase.Moments[moment]) == Moment:
                RHS = RHS - self.LoadCase.Moments[moment].Magnitude
            
            if type(self.LoadCase.Moments[moment]) == ReactionMoment:
                LHS[self.LoadCase.Moments[moment].Magnitude.VariableNumber-1] = 1

        return LHS, RHS

class BeamElement():
    def __init__(self,StartPos,EndPos,LoadCaseEnvironment) -> None:
        self.LoadCase = LoadCaseEnvironment
        self.StartPos = StartPos
        self.EndPos = EndPos
        
        self.ListActiveForces = []
        self.InitializeIValue()
        self.EvaluateDistLoading()

    def InitializeIValue(self):
        self.I_Value = (self.LoadCase.ReturnIValue(self.StartPos)+self.LoadCase.ReturnIValue(self.EndPos))/2

    ### calculating the Distrubuted loading equation for this element
    def EvaluateDistLoading(self):
        Q1 = 0
        Q2 = 0
        x1 = 0
        x2 = 1
        self.GradientQ = 0
        self.Y_InterceptQ = 0
        for dist in self.LoadCase.DistributedLoadings:
            loading = self.LoadCase.DistributedLoadings[dist]
            for idx in range(len(loading.Positions)):
                if idx !=0:
                    if (self.StartPos >= loading.Positions[idx-1]) and (self.StartPos < loading.Positions[idx]):
                        # Linear interpolation
                        Q1 = loading.Magnitudes[idx-1]
                        Q2 = loading.Magnitudes[idx]
                        x1 = loading.Positions[idx-1]
                        x2 = loading.Positions[idx]
                        grad = (Q2-Q1)/(x2-x1)
                        intercept = Q1 - grad*x1
                        self.GradientQ = self.GradientQ + grad
                        self.Y_InterceptQ = self.Y_InterceptQ + intercept

        self.CumulativeShearForce = np.mean([Q2,Q1])*(self.EndPos-self.StartPos)
        if self.CumulativeShearForce !=0:
            self.CGx = (1/self.CumulativeShearForce)*(self.Y_InterceptQ*self.EndPos**2/2 - self.Y_InterceptQ*self.StartPos**2/2 + self.GradientQ*self.EndPos**3/3 - self.GradientQ*self.StartPos**3/3)
        else:
            self.CGx = 0

    ### Still Required for Future implementation;
    # Evaluate Shear force Equation
    # Evaluate Bending moment Equation 
    #                
    def EvaluateAngleEquation(self, Position):
        ### Returns the LHS and RHS of the angle equation.
        RHS = 0
        LHS = np.zeros(len(self.LoadCase.Unknowns))

        
        for i in self.LoadCase.Elements:
            if i == self:             
                RHS = RHS - (self.Y_InterceptQ*Position**3/6 - self.Y_InterceptQ*Position**2*self.StartPos/2 + self.Y_InterceptQ*Position*self.StartPos**2/2 - self.Y_InterceptQ*self.StartPos**3/6 + self.GradientQ*Position**4/24 - self.GradientQ*Position**2*self.StartPos**2/4 + self.GradientQ*Position*self.StartPos**3/3 - self.GradientQ*self.StartPos**4/8)
                break
            else:
                RHS = RHS - (i.Y_InterceptQ*Position**2*i.EndPos/2 - i.Y_InterceptQ*Position**2*i.StartPos/2 - i.Y_InterceptQ*Position*i.EndPos**2/2 + i.Y_InterceptQ*Position*i.StartPos**2/2 + i.Y_InterceptQ*i.EndPos**3/6 - i.Y_InterceptQ*i.StartPos**3/6 + i.GradientQ*Position**2*i.EndPos**2/4 - i.GradientQ*Position**2*i.StartPos**2/4 - i.GradientQ*Position*i.EndPos**3/3 + i.GradientQ*Position*i.StartPos**3/3 + i.GradientQ*i.EndPos**4/8 - i.GradientQ*i.StartPos**4/8)

                
        

        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]

            if type(CurrentForce) == Force:
                if Position>= CurrentForce.Position:
                    RHS = RHS - CurrentForce.Magnitude*(1/2)*(Position-CurrentForce.Position)**2

            if type(CurrentForce) == ReactionForce:
                if Position>= CurrentForce.Position:
                    LHS[CurrentForce.Magnitude.VariableNumber - 1] = (1/2)*(Position-CurrentForce.Position)**2


        for moment in self.LoadCase.Moments:
            CurrentMoment = self.LoadCase.Moments[moment]

            if type(CurrentMoment) == Moment:
                if Position>= CurrentMoment.Position:
                    RHS = RHS + CurrentMoment.Magnitude*(Position-CurrentMoment.Position)

            if type(CurrentMoment) == ReactionMoment:
                if Position>= CurrentMoment.Position:
                    LHS[CurrentMoment.Magnitude.VariableNumber - 1] = -(Position-CurrentMoment.Position)


        LHS[self.LoadCase.C1.VariableNumber - 1] = 1

        return LHS/(self.LoadCase.Lab.Beam.ElasticModulus*self.I_Value), RHS/(self.LoadCase.Lab.Beam.ElasticModulus*self.I_Value)

    def EvaluateDeflectionEquation(self, Position):
        ### Returns the LHS and RHS of the devlection equation.

        RHS = 0
        LHS = np.zeros(len(self.LoadCase.Unknowns))

        
                
        for i in self.LoadCase.Elements:
            if i == self:
                RHS = RHS - (self.Y_InterceptQ*Position**4/24 - self.Y_InterceptQ*Position**3*self.StartPos/6 + self.Y_InterceptQ*Position**2*self.StartPos**2/4 - self.Y_InterceptQ*Position*self.StartPos**3/6 + self.Y_InterceptQ*self.StartPos**4/24 + self.GradientQ*Position**5/120 - self.GradientQ*Position**3*self.StartPos**2/12 + self.GradientQ*Position**2*self.StartPos**3/6 - self.GradientQ*Position*self.StartPos**4/8 + self.GradientQ*self.StartPos**5/30)
                break
            else:
                RHS = RHS - (i.Y_InterceptQ*Position**3*i.EndPos/6 - i.Y_InterceptQ*Position**3*i.StartPos/6 - i.Y_InterceptQ*Position**2*i.EndPos**2/4 + i.Y_InterceptQ*Position**2*i.StartPos**2/4 + i.Y_InterceptQ*Position*i.EndPos**3/6 - i.Y_InterceptQ*Position*i.StartPos**3/6 - i.Y_InterceptQ*i.EndPos**4/24 + i.Y_InterceptQ*i.StartPos**4/24 + i.GradientQ*Position**3*i.EndPos**2/12 - i.GradientQ*Position**3*i.StartPos**2/12 - i.GradientQ*Position**2*i.EndPos**3/6 + i.GradientQ*Position**2*i.StartPos**3/6 + i.GradientQ*Position*i.EndPos**4/8 - i.GradientQ*Position*i.StartPos**4/8 - i.GradientQ*i.EndPos**5/30 + i.GradientQ*i.StartPos**5/30)

        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]

            if type(CurrentForce) == Force:
                if Position>= CurrentForce.Position:
                    RHS = RHS - CurrentForce.Magnitude*(1/6)*(Position-CurrentForce.Position)**3

            if type(CurrentForce) == ReactionForce:
                if Position>= CurrentForce.Position:
                    LHS[CurrentForce.Magnitude.VariableNumber - 1] = (1/6)*(Position-CurrentForce.Position)**3

        for moment in self.LoadCase.Moments:
            CurrentMoment = self.LoadCase.Moments[moment]
            

            if type(CurrentMoment) == Moment:
                if Position>= CurrentMoment.Position:
                    RHS = RHS + CurrentMoment.Magnitude*(1/2)*(Position-CurrentMoment.Position)**2


            if type(CurrentMoment) == ReactionMoment:
                if Position>= CurrentMoment.Position:
                    LHS[CurrentMoment.Magnitude.VariableNumber - 1] = -(1/2)*(Position-CurrentMoment.Position)**2


        LHS[self.LoadCase.C1.VariableNumber - 1] = Position
        LHS[self.LoadCase.C2.VariableNumber - 1] = 1

        return LHS/(self.LoadCase.Lab.Beam.ElasticModulus*self.I_Value), RHS/(self.LoadCase.Lab.Beam.ElasticModulus*self.I_Value)

    ### Plotting Functions:

    def DistributedLoadingFunction(self,xs):
        return self.GradientQ*xs + self.Y_InterceptQ

    def ShearForceFunction(self,xs):
        ForcePresent = np.zeros(len(xs))
        ## Distributed loading:

        
                
        for i in self.LoadCase.Elements:
            if i == self:
                ForcePresent = ForcePresent + (self.Y_InterceptQ*xs - self.Y_InterceptQ*self.StartPos + self.GradientQ*xs**2/2 - self.GradientQ*self.StartPos**2/2)
                break
            else:
                ForcePresent = ForcePresent + (i.Y_InterceptQ*i.EndPos - i.Y_InterceptQ*i.StartPos + i.GradientQ*i.EndPos**2/2 - i.GradientQ*i.StartPos**2/2)


        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]
            if CurrentForce.Position <= self.StartPos:
                if type(CurrentForce) == Force:
                    ForcePresent = ForcePresent + CurrentForce.Magnitude
                if type(CurrentForce) == ReactionForce:
                    ForcePresent = ForcePresent + CurrentForce.Magnitude.Magnitude

        return ForcePresent

    def BendingMomentFunction(self,xs):
        MomentPresent = np.zeros(len(xs))
        ## Distributed loading:

        
                
        for i in self.LoadCase.Elements:
            if i == self:
                MomentPresent = MomentPresent + (self.Y_InterceptQ*xs**2/2 - self.Y_InterceptQ*xs*self.StartPos + self.Y_InterceptQ*self.StartPos**2/2 + self.GradientQ*xs**3/6 - self.GradientQ*xs*self.StartPos**2/2 + self.GradientQ*self.StartPos**3/3)
                break
            else:
                # print(i.CGx)
                MomentPresent = MomentPresent + (i.Y_InterceptQ*xs*i.EndPos - i.Y_InterceptQ*xs*i.StartPos - i.Y_InterceptQ*i.EndPos**2/2 + i.Y_InterceptQ*i.StartPos**2/2 + i.GradientQ*xs*i.EndPos**2/2 - i.GradientQ*xs*i.StartPos**2/2 - i.GradientQ*i.EndPos**3/3 + i.GradientQ*i.StartPos**3/3)

        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]
            if CurrentForce.Position <= self.StartPos:
                if type(CurrentForce) == Force:
                    MomentPresent = MomentPresent + (CurrentForce.Magnitude*(xs-CurrentForce.Position))
                if type(CurrentForce) == ReactionForce:
                    MomentPresent = MomentPresent + (CurrentForce.Magnitude.Magnitude*(xs-CurrentForce.Position))

        for moment in self.LoadCase.Moments:
            CurrentMoment = self.LoadCase.Moments[moment]
            if CurrentMoment.Position <= self.StartPos:
                if type(CurrentMoment) == Moment:
                    MomentPresent = MomentPresent - CurrentMoment.Magnitude
                if type(CurrentMoment) == ReactionMoment:
                    MomentPresent = MomentPresent - CurrentMoment.Magnitude.Magnitude

        return MomentPresent

    def AngleThetaFunction(self,xs):
        EIs = (self.LoadCase.Lab.Beam.ReturnIvalue(xs)*self.LoadCase.Lab.Beam.ElasticModulus)
        ThetaPresent = np.zeros(len(xs)) 
        ## Distributed loading:

        for i in self.LoadCase.Elements:
            if i == self:
                ThetaPresent = ThetaPresent + (self.Y_InterceptQ*xs**3/6 - self.Y_InterceptQ*xs**2*self.StartPos/2 + self.Y_InterceptQ*xs*self.StartPos**2/2 - self.Y_InterceptQ*self.StartPos**3/6 + self.GradientQ*xs**4/24 - self.GradientQ*xs**2*self.StartPos**2/4 + self.GradientQ*xs*self.StartPos**3/3 - self.GradientQ*self.StartPos**4/8)
                break
            else:
                ThetaPresent = ThetaPresent + (i.Y_InterceptQ*xs**2*i.EndPos/2 - i.Y_InterceptQ*xs**2*i.StartPos/2 - i.Y_InterceptQ*xs*i.EndPos**2/2 + i.Y_InterceptQ*xs*i.StartPos**2/2 + i.Y_InterceptQ*i.EndPos**3/6 - i.Y_InterceptQ*i.StartPos**3/6 + i.GradientQ*xs**2*i.EndPos**2/4 - i.GradientQ*xs**2*i.StartPos**2/4 - i.GradientQ*xs*i.EndPos**3/3 + i.GradientQ*xs*i.StartPos**3/3 + i.GradientQ*i.EndPos**4/8 - i.GradientQ*i.StartPos**4/8)

        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]
            if CurrentForce.Position <= self.StartPos:
                if type(CurrentForce) == Force:
                    ThetaPresent = ThetaPresent + ((CurrentForce.Magnitude/2)*(xs-CurrentForce.Position)**2)    
                if type(CurrentForce) == ReactionForce:
                    ThetaPresent = ThetaPresent + ((CurrentForce.Magnitude.Magnitude/2)*(xs-CurrentForce.Position)**2)

        for moment in self.LoadCase.Moments:
            CurrentMoment = self.LoadCase.Moments[moment]
            if CurrentMoment.Position <= self.StartPos:
                if type(CurrentMoment) == Moment:
                    ThetaPresent = ThetaPresent - (CurrentMoment.Magnitude*(xs-CurrentMoment.Position))
                if type(CurrentMoment) == ReactionMoment:
                    ThetaPresent = ThetaPresent - (CurrentMoment.Magnitude.Magnitude*(xs-CurrentMoment.Position))


        return (ThetaPresent + (self.LoadCase.C1.Magnitude))*(180/np.pi)/EIs

    def DisplacementFunction(self,xs):
        EIs = (self.LoadCase.Lab.Beam.ReturnIvalue(xs)*self.LoadCase.Lab.Beam.ElasticModulus)
        DisplacePresent = np.zeros(len(xs))
        ## Distributed loading:

        for i in self.LoadCase.Elements:
            if i == self:
                DisplacePresent = DisplacePresent + (self.Y_InterceptQ*xs**4/24 - self.Y_InterceptQ*xs**3*self.StartPos/6 + self.Y_InterceptQ*xs**2*self.StartPos**2/4 - self.Y_InterceptQ*xs*self.StartPos**3/6 + self.Y_InterceptQ*self.StartPos**4/24 + self.GradientQ*xs**5/120 - self.GradientQ*xs**3*self.StartPos**2/12 + self.GradientQ*xs**2*self.StartPos**3/6 - self.GradientQ*xs*self.StartPos**4/8 + self.GradientQ*self.StartPos**5/30)
                break
            else:
                DisplacePresent = DisplacePresent + (i.Y_InterceptQ*xs**3*i.EndPos/6 - i.Y_InterceptQ*xs**3*i.StartPos/6 - i.Y_InterceptQ*xs**2*i.EndPos**2/4 + i.Y_InterceptQ*xs**2*i.StartPos**2/4 + i.Y_InterceptQ*xs*i.EndPos**3/6 - i.Y_InterceptQ*xs*i.StartPos**3/6 - i.Y_InterceptQ*i.EndPos**4/24 + i.Y_InterceptQ*i.StartPos**4/24 + i.GradientQ*xs**3*i.EndPos**2/12 - i.GradientQ*xs**3*i.StartPos**2/12 - i.GradientQ*xs**2*i.EndPos**3/6 + i.GradientQ*xs**2*i.StartPos**3/6 + i.GradientQ*xs*i.EndPos**4/8 - i.GradientQ*xs*i.StartPos**4/8 - i.GradientQ*i.EndPos**5/30 + i.GradientQ*i.StartPos**5/30)

        for force in self.LoadCase.Forces:
            CurrentForce = self.LoadCase.Forces[force]
            if CurrentForce.Position <= self.StartPos:
                if type(CurrentForce) == Force:
                    DisplacePresent = DisplacePresent + ((CurrentForce.Magnitude/6)*(xs-CurrentForce.Position)**3)  
                if type(CurrentForce) == ReactionForce:
                    DisplacePresent = DisplacePresent + ((CurrentForce.Magnitude.Magnitude/6)*(xs-CurrentForce.Position)**3)

        for moment in self.LoadCase.Moments:
            CurrentMoment = self.LoadCase.Moments[moment]
            if CurrentMoment.Position <= self.StartPos:
                if type(CurrentMoment) == Moment:
                    DisplacePresent = DisplacePresent - ((CurrentMoment.Magnitude/2)*(xs-CurrentMoment.Position)**2)
                if type(CurrentMoment) == ReactionMoment:
                    DisplacePresent = DisplacePresent - ((CurrentMoment.Magnitude.Magnitude/2)*(xs-CurrentMoment.Position)**2)


        print
        return (DisplacePresent + (xs*self.LoadCase.C1.Magnitude) + (self.LoadCase.C2.Magnitude))/EIs

class UnknownConstant():
    def __init__(self,Name,LoadCase) -> None:
        self.Name = Name
        self.Magnitude = 0
        self.LoadCase = LoadCase
        self.LoadCase.Unknowns.append(self)
        self.VariableNumber = len(self.LoadCase.Unknowns)

