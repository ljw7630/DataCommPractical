import receiver_preprocess
import os

num_of_pairs = range(8, 26, 2)
num_of_times = 5

def main():
	for num_of_pair in num_of_pairs:
		for i in range(0, num_of_times):
			os.system("run.bat " + str(num_of_pair) + " 0")

if __name__ == '__main__':
	main()