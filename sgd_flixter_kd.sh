export GRAPHCHI_ROOT=../
for m in '60' '70' '80' '90' '100'
do
	for l in  '1e-3'
  	do
		../bin/myapps/sgd_static_kdtree --D=$m --training=../../data/flixter/ratings_flixter_mm --result=../../result/gc_flixter_res_kd \
		  --sgd_gamma=1e-4 --sgd_lambda=1e-4 --max_iter=100 --n_top=50 --quiet=1 --tol=1e-3
    	done
done
