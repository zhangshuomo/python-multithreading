# synchronization of n threads 
# i.e. guarantee the exec times of all the n threads are the same
import threading as td
from threading import Condition
import time
import numpy as np

n_of_thread = 100
loop_counter=dict(zip(['Thread-{}'.format(i) for i in range(1, n_of_thread+1)], [0,]*n_of_thread))	# a counter for each thread exec times
cond = Condition()

def func():
	while True:
		cond.acquire()
		values = list(loop_counter.values())
		called_times = loop_counter[td.current_thread().name]
		if called_times > min(values):	cond.wait()	# there exist other threads whose exec times is smaller than the current thread
		loop_counter[td.current_thread().name] += 1
		time.sleep(abs(np.random.randn())*0.00001)
		if len(set(list(loop_counter.values())))==1: # all elements in call_counter dict are the same, then notify all the waiting thread
			print(loop_counter)
			cond.notifyAll()
		cond.release()
		
for i in range(n_of_thread):
	t = td.Thread(target=func)
	t.start()
