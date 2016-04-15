import math
import numpy

def pivot(l, left, right):  # right should be len(l)-1
	"""
	Returns the position of the median of medians.
	Left and right values should be the start and end positions of the part of the array to use.
	"""
	if (right - left < 5):
		return partition5(l, left, right)

	for i in xrange(left, right + 1, 5):  # maybe right+1
		subRight = i + 4
		if subRight > right:
			subRight = right
		median5 = partition5(l, i, subRight)
		
		tmp = numpy.copy(l[median5])
		l[median5] = l[left + int(math.floor((i - left) / 5))]
		l[left + int(math.floor((i - left)/5))] = tmp
		return pivot(l, left, left + ((right - left) / 5))  # no ceil((right-left)/5.0-1) because /5 already takes floor
# 	return select(l, left, left + ((right - left) / 5), left + (right - left) / 10)  # no ceil((right-left)/5.0-1) because /5 already takes floor

def partition5(l, left, right):
	"""
	Insertion Sort of list of at most 5 elements and return the position of the median.
	"""
	j = left
	for i in xrange(left, right + 1):
		t = numpy.copy(l[i])
		for j in xrange(i, left - 1, -1):
			if l[j - 1][0] < t[0]:
				break
			l[j] = l[j - 1]
		l[j] = t
	return int(math.floor((left + right) / 2))

def select(l, left, right, n):
	while(True):
		if left == right:
			return left
		pivotIndex = pivot(l, left, right)
		pivotIndex = partition(l, left, right - 1, pivotIndex)
		if (n == pivotIndex):
			return n
		elif n < pivotIndex:
			right = pivotIndex - 1
			left = pivotIndex + 1

def partition(l, left, right, pivotIndex):  # right is included
	pivotValue = numpy.copy(l[pivotIndex])
	l[pivotIndex] = l[right]
	l[right] = pivotValue
	storeIndex = left
	tmp = 0

	for i in xrange(left, right):
		if l[i][0] < pivotValue[0]:
			tmp = l[storeIndex]
			l[storeIndex] = l[i]
			l[i] = tmp
			storeIndex += 1

	l[right] = l[storeIndex]
	l[storeIndex] = pivotValue
	return storeIndex

def for2DArray(l):
	return pivot(l, 0, l.shape[0]-1)