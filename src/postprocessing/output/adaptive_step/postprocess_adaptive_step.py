import numpy as np 
import matplotlib.pyplot as plt

def read_stats_data(problem, methods, tols_dict, controller):
	data = {}
	for method in methods:
		data[method] = {}
		for tol in tols_dict[method]:
			#print("\tMethod: " + method, flush=True)
			stats_filename = "./output/" + problem + "/" + problem + "_AdaptiveStep_" + tol + "_" + controller + "_LASA-mean_" + method + "_stats.csv"
			
			#print(stats_filename)
			try:
				stats_data = np.genfromtxt(stats_filename, delimiter=",")
				#print(stats_data)

				method_data = { 
					"total_timesteps": stats_data[0],
					"total_successful_timesteps": stats_data[1],
					"total_microtimesteps": stats_data[2],
					"total_successful_microtimesteps": stats_data[3],
					"rel_err": stats_data[4],
					"abs_err": stats_data[5],
					"fast_function_evals": stats_data[7],
					"slow_function_evals": stats_data[8],
					"implicit_function_evals": stats_data[9],
					"explicit_function_evals": stats_data[10],
					"fast_jacobian_evals": stats_data[12],
					"slow_jacobian_evals": stats_data[13],
					"implicit_jacobian_evals": stats_data[14],
					"slow_nonlinear_solves": stats_data[15],
					"runtime": stats_data[16],
					"status": stats_data[17]
				};

				data[method][tol] = method_data
			except:
				print("Problem: (" + problem + ") has no (stats) data from method: (" + method + ")", flush=True)
	return data

def plot_fast_evals(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["x"].append(stats_data_dict[method][tol]["fast_function_evals"])
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.loglog(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=3, borderaxespad=0.)
	#ax1.legend(methods,fontsize=14)
	ax1.set_ylabel("Max Error",fontsize=14)
	ax1.set_xlabel("Fast Function Calls",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_fastevals.png")
	
def plot_slow_evals(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["x"].append(stats_data_dict[method][tol]["slow_function_evals"])
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.semilogx(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	#ax1.legend(methods,fontsize=14)
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
		   ncol=3, borderaxespad=0.)
	ax1.set_ylabel("Max Error",fontsize=14)
	ax1.set_xlabel("Slow Function Calls",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_slowevals.png")
	
def plot_implicit_evals(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["x"].append(stats_data_dict[method][tol]["implicit_function_evals"])
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.loglog(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	#ax1.legend(methods,fontsize=14)
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
		   ncol=3, borderaxespad=0.)
	ax1.set_ylabel("Max Error",fontsize=14)
	ax1.set_xlabel("Implicit Function Calls",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_implicitevals.png")
	
def plot_implicit_solves(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["x"].append(stats_data_dict[method][tol]["slow_nonlinear_solves"])
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.loglog(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	#ax1.legend(methods,fontsize=14)
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
		   ncol=3, borderaxespad=0.)
	ax1.set_ylabel("Max Error",fontsize=14)
	ax1.set_xlabel("Implicit Solves",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_implicitsolves.png")
	
def plot_err_deviation(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["x"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"]/float(tol))
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.semilogx(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	#ax1.legend(methods,fontsize=14)
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
		   ncol=3, borderaxespad=0.)
	ax1.set_xlabel("Max Error",fontsize=14)
	ax1.set_ylabel("Error Deviation",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_errdev.png")

def plot_err(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["x"].append(float(tol))
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
		
	print(plot_data)
	fig,ax1 = plt.subplots()
	for method in methods:
		ax1.semilogx(plot_data[method]["x"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"])
	#ax1.legend(methods,fontsize=14)
	ax1.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
		   ncol=3, borderaxespad=0.)
	ax1.set_xlabel("Max Error",fontsize=14)
	ax1.set_ylabel("Error Deviation",fontsize=14)
	ax1.tick_params(labelsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_" + method + "_err.png")

def plot_two_factors(problem, methods, controller, tols_dict, stats_data_dict, styles_dict):
	plot_data = {}
	legend_items = []
	for method in methods:
		plot_data[method] = {"x1": [], "x2": [], "y": []}
		legend_items.append(styles_dict[method]["legend_label"])
		for tol in tols_dict[method]:
			plot_data[method]["y"].append(stats_data_dict[method][tol]["abs_err"])
			plot_data[method]["x1"].append(stats_data_dict[method][tol]["fast_function_evals"])
			plot_data[method]["x2"].append(stats_data_dict[method][tol]["slow_nonlinear_solves"])
		
	print(plot_data)
	fig,(ax1,ax2) = plt.subplots(1,2,figsize=(7.6,5.7))
	for method in methods:
		ax1.loglog(plot_data[method]["x1"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"],markersize=7)
		ax2.loglog(plot_data[method]["x2"],plot_data[method]["y"],linestyle='-',marker=styles_dict[method]["marker"], color=styles_dict[method]["color"],markersize=7)
	
	fig.legend(legend_items,prop={'size':11},bbox_to_anchor=(0.15,0.77,0.815,0.2), mode="expand", loc="lower left",
               ncol=3, borderaxespad=0.)
	ax1.set_box_aspect(1)
	ax2.set_box_aspect(1)
	ax1.set_ylabel("Max Error",fontsize=11)
	ax1.set_xlabel("Fast Function Evaluations",fontsize=11)
	ax2.set_ylabel("Max Error",fontsize=11)
	ax2.set_xlabel("Implicit Solves",fontsize=11)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/adaptive_step/plots/" + problem + "_twofactors.pdf",bbox_inches="tight")

def main():
	problem = "Brusselator1DIMEX"
	methods = ["IMEXMRISR21", "IMEXMRISR32"]
	tols_dict = {
		"IMEXMRISR21": ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9"],
		"IMEXMRISR32": ["1e-1", "1e-2", "1e-4", "1e-6", "1e-7", "1e-8", "1e-9"]
	}
	styles_dict = {
		"IMEXMRISR21": {
			"color": "blue",
			"marker": '^',
			"legend_label": "IMEX-MRI-SR2(1)"
		},
		"IMEXMRISR32": {
			"color": "orange",
			"marker": 'o',
			"legend_label": "IMEX-MRI-SR3(2)"
		}
	}
	controller = "ConstantConstant"

	stats_data = read_stats_data(problem, methods, tols_dict, controller)
	
	print(stats_data)

	plot_fast_evals(problem, methods, controller, tols_dict, stats_data, styles_dict)
	plot_slow_evals(problem, methods, controller, tols_dict, stats_data, styles_dict)
	plot_implicit_evals(problem, methods, controller, tols_dict, stats_data, styles_dict)
	plot_implicit_solves(problem, methods, controller, tols_dict, stats_data, styles_dict)
	plot_err_deviation(problem, methods, controller, tols_dict, stats_data, styles_dict)
	plot_two_factors(problem, methods, controller, tols_dict, stats_data, styles_dict)


if __name__ == "__main__":
	main()