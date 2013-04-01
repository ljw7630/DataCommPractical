import os

result_files = ["rateAverageSum1Result.txt" \
	, "rateAverageVariance1Result.txt" \
	, "rateIndexSum1Result.txt" \
	, "rateWeightSum1Result.txt"]

def main():
	all_files = get_experiment_files()

	for to_be_processed in files:
		files = [f for f in all_files if f.endswith(to_be_processed)]
		print to_be_processed, files

def get_experiment_files():
	first_level_dirs = get_immediate_subdirectories('.')
	print first_level_dirs
	for first_level_dir in first_level_dirs:
		second_level_dirs = get_immediate_subdirectories(first_level_dir)
		print "second level", second_level_dirs
		for second_level_dir in second_level_dirs:
			files = get_immediate_files(second_level_dir)
	return files

def get_immediate_files(path):
	return [os.path.join(path, name) for name in os.listdir(path) if os.path.isfile(os.path.join(path, name))]

def get_immediate_subdirectories(path):
	return [os.path.join(path, name) for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]

if __name__ == '__main__':
	main()
