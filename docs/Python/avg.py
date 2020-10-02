def average(list):
    count=0
    sum  = 0.0
    for item in list:
    	count += 1
	sum   += item
    return sum / count

print(average([100,200,300,150,110,99]))
