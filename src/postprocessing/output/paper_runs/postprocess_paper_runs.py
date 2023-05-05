import numpy as np 
import matplotlib.pyplot as plt

def read_stats_data(problem, methods, m):
	data = {}
	for method in methods:
		stats_filename = "./output/" + problem + "/" + problem + "_PaperRuns_" + method + "_M" + str(m) + "_stats.csv"
		try:
			stats_data = np.genfromtxt(stats_filename, delimiter=",")
			method_data = { 
				"hs": stats_data[:,0],
				"errs": stats_data[:,1],
				"fast_function_evals": stats_data[:,2],
				"slow_function_evals": stats_data[:,3],
				"implicit_function_evals": stats_data[:,4],
				"explicit_function_evals": stats_data[:,5],
				"fast_jacobian_evals": stats_data[:,6],
				"slow_jacobian_evals": stats_data[:,7],
				"implicit_jacobian_evals": stats_data[:,8],
				"slow_nonlinear_solves": stats_data[:,9],
				"runtimes": stats_data[:,10],
				"m": m
			};
			data[method] = method_data
		except:
			print("Problem: (" + problem + ") has no data from method: (" + method + ") for m: (" + str(m) + ")\n")
	return data

def main():
	problems = ['KPR']
	imex_methods = ['LieTrotter', 'StrangMarchuk', 'IMEXMRISR21', 'IMEXMRISR32', 'IMEXMRI3a', 'IMEXMRI3b']
	imex_methods_pretty = ['Lie-Trotter','Strang-Marchuk','IMEX-MRI-SR2(1)','IMEX-MRI-SR3(2)','IMEX-MRI-GARK3a', 'IMEX-MRI-GARK3b']
	m = 10
	imex_markers = ['<','^','>','v','d','s','.','o']
	imex_colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'k', 'gray']
	marker_size = 7
	
	imex_data = read_stats_data("KPR", imex_methods, m)
		
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(imex_methods)):
		method = imex_methods[i]
		method_data = imex_data[method]
		plt.loglog(method_data["hs"],method_data["errs"],marker=imex_markers[i],color=imex_colors[i],markersize=marker_size)
		p = np.polyfit(np.log10(method_data["hs"]),np.log10(method_data["errs"]),1)
		legend_items.append(imex_methods_pretty[i] + " (" + str(np.round(p[0],2)) + ")")
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.gca().invert_xaxis()
	plt.xlabel("H",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	#plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/KPR_convergence_imex.pdf", bbox_inches="tight")
	plt.close()
			
	
	merk_methods = ['IMEXMRISRMERK2', 'IMEXMRISRMERK3', 'IMEXMRISRMERK4', 'IMEXMRISRMERK5']
	merk_methods_pretty = ['MERK2', 'MERK3', 'MERK4', 'MERK5']
	m = 10
	merk_markers = ['<','>','^','v']
	merk_colors = ['blue', 'orange', 'green', 'red']

	merk_data = read_stats_data('KPR', merk_methods, m)
		
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(merk_methods)):
		method = merk_methods[i]
		method_data = merk_data[method]
		plt.loglog(method_data["hs"],method_data["errs"],marker=merk_markers[i],color=merk_colors[i],markersize=marker_size)
		p = np.polyfit(np.log10(method_data["hs"]),np.log10(method_data["errs"]),1)
		legend_items.append(merk_methods_pretty[i] + " (" + str(np.round(p[0],2)) + ")")
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.gca().invert_xaxis()
	plt.xlabel("H",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	#plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/KPR_convergence_merk.pdf", bbox_inches="tight")
	plt.close()
	
	
	im_methods = ['LieTrotter', 'StrangMarchuk', 'IMEXMRISR21', 'IMEXMRISR32', 'IMEXMRI3a', 'IMEXMRI3b']
	im_methods_pretty = ['Lie-Trotter','Strang-Marchuk','IMEX-MRI-SR2(1)','IMEX-MRI-SR3(2)','IMEX-MRI-GARK3a','IMEX-MRI-GARK3b', 'MRI-GARK21a', 'MRI-GARK34a']
	m = 10
	im_markers = ['<','^','>','v','d','s','.','o','p','P','*']
	im_data_201 = read_stats_data('BrusselatorPDE201', im_methods, m)
		
	im_colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'k', 'gray', 'pink', 'y', 'cyan']
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(im_methods)):
		method = im_methods[i]
		method_data = im_data_201[method]
		plt.loglog(method_data["runtimes"]/1000,method_data["errs"],marker=im_markers[i],color=im_colors[i],markersize=marker_size)
		legend_items.append(im_methods_pretty[i])
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.xlabel("Runtime (s)",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE201_runtimes_im.pdf")
	plt.close()
	
	im_data_801 = read_stats_data('BrusselatorPDE801', im_methods, m)
		
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(im_methods)):
		method = im_methods[i]
		method_data = im_data_801[method]
		plt.loglog(method_data["runtimes"]/1000,method_data["errs"],marker=im_markers[i],color=im_colors[i],markersize=marker_size)
		legend_items.append(im_methods_pretty[i])
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.xlabel("Runtime (s)",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE801_runtimes_im.pdf")
	plt.close()
	
	fig, (ax1, ax2) = plt.subplots(1,2)
	legend_items = []
	for i in range(len(im_methods)):
		method = im_methods[i]
		method_data_201 = im_data_201[method]
		method_data_801 = im_data_801[method]
		ax1.loglog(method_data_201["runtimes"]/1000,method_data_201["errs"],marker=im_markers[i],color=im_colors[i],markersize=marker_size)
		ax2.loglog(method_data_801["runtimes"]/1000,method_data_801["errs"],marker=im_markers[i],color=im_colors[i],markersize=marker_size)
		legend_items.append(im_methods[i])
	
	ax1.text(0.925,0.925,'201 grid',ha='right',va='top',transform=ax1.transAxes,bbox={'facecolor':'None','alpha':0.5,'pad':5})
	ax2.text(0.925,0.925,'801 grid',ha='right',va='top',transform=ax2.transAxes,bbox={'facecolor':'None','alpha':0.5,'pad':5})
	ax1.set_box_aspect(1)
	ax2.set_box_aspect(1)
	ax1.set_ylabel("Max Error",fontsize=14)
	ax2.set_ylabel("Max Error",fontsize=14)
	ax1.set_xlabel("Runtime (s)",fontsize=14)
	ax2.set_xlabel("Runtime (s)",fontsize=14)
	#plt.tick_params(axis='both',which='major',labelsize=12)
	fig.legend(legend_items,prop={'size':10},bbox_to_anchor=(0.15,0.75,0.815,0.2), mode="expand", loc="lower left",
               ncol=3, borderaxespad=0.)
	plt.tight_layout()
	#ax1.set_xticks(fontsize=12)
	#ax1.set_yticks(fontsize=12)
	#ax2.set_xticks(fontsize=12)
	#ax2.set_yticks(fontsize=12)
	
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE201801_runtimes_im.pdf",bbox_inches="tight")
	plt.close()
	
	
	imex_methods = ['LieTrotter', 'StrangMarchuk', 'IMEXMRISR21', 'IMEXMRISR32', 'IMEXMRI3a', 'IMEXMRI3b']
	imex_methods_pretty = ['Lie-Trotter','Strang-Marchuk','IMEX-MRI-SR2(1)','IMEX-MRI-SR3(2)','IMEX-MRI-GARK3a', 'IMEX-MRI-GARK3b']
	m = 10
	imex_markers = ['<','^','>','v','d','s','.','o']
	imex_data_201 = read_stats_data('BrusselatorPDE201', imex_methods, m)
		
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(imex_methods)):
		method = imex_methods[i]
		method_data = imex_data_201[method]
		plt.loglog(method_data["runtimes"]/1000,method_data["errs"],marker=imex_markers[i],markersize=marker_size)
		legend_items.append(im_methods_pretty[i])
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.xlabel("Runtime (s)",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE201_runtimes_imex.pdf")
	plt.close()
	
	imex_data_801 = read_stats_data('BrusselatorPDE801', imex_methods, m)
		
	fig = plt.figure()
	#plt.title("H vs. Err for " + problem + " Problem with M=" + str(m))
	legend_items = []
	for i in range(len(imex_methods)):
		method = imex_methods[i]
		method_data = imex_data_801[method]
		plt.loglog(method_data["runtimes"]/1000,method_data["errs"],marker=imex_markers[i],markersize=marker_size)
		legend_items.append(im_methods_pretty[i])
	plt.legend(legend_items,prop={'size':11},bbox_to_anchor=(0,1.02,1,0.2), mode="expand", loc="lower left",
               ncol=2, borderaxespad=0.0)
	plt.xlabel("Runtime (s)",fontsize=14)
	plt.ylabel("Max Error",fontsize=14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.tight_layout()
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE801_runtimes_imex.pdf")
	plt.close()
	
	imex_colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'k', 'gray']
	fig, (ax1, ax2) = plt.subplots(1,2,figsize=(7.6,5.7))
	#fig = plt.figure(figsize=(7.2,5.4))
	#ax1 = fig.sub
	legend_items = []
	for i in range(len(imex_methods)):
		method = imex_methods[i]
		method_data_201 = im_data_201[method]
		method_data_801 = im_data_801[method]
		ax1.loglog(method_data_201["runtimes"]/1000,method_data_201["errs"],marker=imex_markers[i],color=imex_colors[i],markersize=marker_size)
		ax2.loglog(method_data_801["runtimes"]/1000,method_data_801["errs"],marker=imex_markers[i],color=imex_colors[i],markersize=marker_size)
		legend_items.append(imex_methods_pretty[i])
	
	#ax1.legend(labels=None,title='201 grid')
	#ax2.legend(labels=None,title='801 grid')
	ax1.text(0.925,0.925,'201 grid',ha='right',va='top',transform=ax1.transAxes,bbox={'facecolor':'None','alpha':0.5,'pad':5})
	ax2.text(0.925,0.925,'801 grid',ha='right',va='top',transform=ax2.transAxes,bbox={'facecolor':'None','alpha':0.5,'pad':5})
	ax1.set_box_aspect(1)
	ax2.set_box_aspect(1)
	ax1.set_ylabel("Max Error",fontsize=11)
	ax2.set_ylabel("Max Error",fontsize=11)
	ax1.set_xlabel("Runtime (s)",fontsize=11)
	ax2.set_xlabel("Runtime (s)",fontsize=11)
	#plt.tick_params(axis='both',which='major',labelsize=12)
	fig.legend(legend_items,prop={'size':11},bbox_to_anchor=(0.15,0.77,0.815,0.2), mode="expand", loc="lower left",
               ncol=3, borderaxespad=0.)
	plt.tight_layout()
	#ax1.set_xticks(fontsize=12)
	#ax1.set_yticks(fontsize=12)
	#ax2.set_xticks(fontsize=12)
	#ax2.set_yticks(fontsize=12)
	
	fig.savefig("./postprocessing/output/paper_runs/plots/BrusselatorPDE201801_runtimes_imex.pdf",bbox_inches="tight")
	plt.close()
	


if __name__ == "__main__":
	main()