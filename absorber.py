import pyprop.core as core
from numpy import ones, double, cos, pi

from utils import RegisterProjectNamespace, RegisterAll

@RegisterProjectNamespace
@RegisterAll
class GlobalCosineAbsorber(core.AbsorberModel):

	def ApplyConfigSection(self, config):
		self.Strength = config.strength

	def SetupStep(self, localGrid, globalGrid):
		#get leftmost global grid position
		minPos = min(globalGrid)
		maxPos = max(globalGrid)

		#get absorber strength
		strength = self.Strength

		#calculate width of grid
		width = maxPos - minPos

		self.Scaling = ones(len(localGrid), dtype=double)

		for i, pos in enumerate(localGrid):	
			self.Scaling[i] = 1.0 - cos(pi * (pos - minPos) / width)**strength

	def	GetScaling(self):
		return self.Scaling


