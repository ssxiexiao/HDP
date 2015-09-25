# -*-coding:utf-8-*-

def sampling(p):
	for i in range(1, len(p)):
		p[i] += p[i-1]
	u = random.random()
	for i in range(len(p)):
		if p[i] > u:
			return i

#Chinese restaurant franchise
class HDP:
	def __init__(self):
		self.group_num = 1
		self.gama = []
		self.alpha = []
		self.data = []
		self.alpha_a = 2
		self.alpha_b = 1
		self.gama_a = 2
		self.gama_b = 1
		# n[j][t] denotes the number of customers sitting at table t in restaurant j
		self.n = []
		# m[j][k] denotes the number of tables in restaurant j serving dish k
		self.m = []
		# k[j][t] denotes the category of dish served by table t in restaurant j 
		self.k = []
		# global_k denotes the menu of all the restaurants
		self.global_k = []

		self.table = []
		self.z = []

	def group_num_setting(self, num):
		self.group_num = num

	def marginalX(self, j, i, specific_k):	

	def pTnew(self, j, i):

	def samplingT(self, j, i):
		p = []
		for t in range(len(self.n[j])):
			result = self.n[j][t] * self.marginalX(j, i, self.k[j][t])
			p.append(result)
		result = self.alpha[j] * self.pTnew(j, i)	
		p.append(result)
		resultT = sampling(p)
		return resultT

	def samplingZ(self, j, i):

	def samplingAlpha(self, j, i):

	def samplingGama(self, j, i):
	
	def _delete(self, j, i):
		z_old = self.z[j][i]
		t_old = self.table[j][i]
		self.n[j][t_old] -= 1


