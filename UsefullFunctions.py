

    def BuildEqSumOfForces(self):
        RHS = 0
        LHS = np.zeros(self.NumOfUnknowns)
        for dist in self.DistributedLoadings:
            for posdx, pos in enumerate(self.DistributedLoadings[dist].Positions):
                if posdx!=0:
                    RHS = RHS - np.mean([self.DistributedLoadings[dist].Magnitudes[posdx-1],self.DistributedLoadings[dist].Magnitudes[posdx]])*(self.DistributedLoadings[dist].Positions[posdx]-self.DistributedLoadings[dist].Positions[posdx-1])

        for force in self.Forces:
            if type(self.Forces[force]) == Force:
                RHS = RHS - self.Forces[force].Magnitude
            
            if type(self.Forces[force]) == ReactionForce:
                LHS[self.Forces[force].VariableNumber-1] = 1

        self.Matrix[0] = LHS
        self.Vector[0] = RHS

    def BuildEqSumOfMoments(self):
        RHS = 0
        LHS = np.zeros(self.NumOfUnknowns)
        for dist in self.DistributedLoadings:
            for posdx, pos in enumerate(self.DistributedLoadings[dist].Positions):
                if posdx!=0:
                    
                    Moment1 = self.DistributedLoadings[dist].Magnitudes[posdx-1]*(self.DistributedLoadings[dist].Positions[posdx]-self.DistributedLoadings[dist].Positions[posdx-1])*np.mean([self.DistributedLoadings[dist].Positions[posdx],self.DistributedLoadings[dist].Positions[posdx-1]])
                    Moment2 = 0.5*(self.DistributedLoadings[dist].Magnitudes[posdx]-self.DistributedLoadings[dist].Magnitudes[posdx-1])*(self.DistributedLoadings[dist].Positions[posdx]-self.DistributedLoadings[dist].Positions[posdx-1])*(self.DistributedLoadings[dist].Positions[posdx-1] + (2/3)*(self.DistributedLoadings[dist].Positions[posdx]-self.DistributedLoadings[dist].Positions[posdx-1]))
                    RHS = RHS - Moment1 - Moment2

        for force in self.Forces:
            if type(self.Forces[force]) == Force:
                RHS = RHS - self.Forces[force].Magnitude*self.Forces[force].Position
            
            if type(self.Forces[force]) == ReactionForce:
                LHS[self.Forces[force].VariableNumber-1] = self.Forces[force].Position

        for moment in self.Moments:
            if type(self.Moments[moment]) == Moment:
                RHS = RHS - self.Moments[moment].Magnitude
            
            if type(self.Moments[moment]) == ReactionMoment:
                LHS[self.Moments[moment].VariableNumber-1] = 1

        self.Matrix[1] = LHS
        self.Vector[1] = RHS