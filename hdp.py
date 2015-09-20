# -*-coding:utf-8-*-

#Chinese restaurant franchise
class HDP:
	def __init__(self):
		self.group_num = 1
		self.gama = []
		self.alpha = []
		# n[j][t] denotes the number of customers sitting at table t in restaurant j
		self.n = []
		# m[j][k] denotes the number of tables in restaurant j serving dish k
		self.m = []
		# k[j][t] denotes the category of dish searved by table t in restaurant j 
		self.k = []
		# global_k denotes the menu of all the restaurants
		self.global_k = []

	def group_num_setting(self, num):
		self.group_num = num

	def samplingT(self, j, i):


