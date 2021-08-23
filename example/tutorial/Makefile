CXX  = g++

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += -std=c++14
CXXFLAGS += -O3

BOOST_PYTHON_LIB = /home/burmist/home2/training/032_boost/boost_1_76_0_install/lib

PYTHON_INCLUDE = /home/burmist/miniconda3/envs/boost_python/include/python3.7m
BOOST_PYTHON_INCLUDE = /home/burmist/home2/training/032_boost/boost_1_76_0_install/include

INCLUDE = -I$(PYTHON_INCLUDE) -I$(BOOST_PYTHON_INCLUDE)
LIB = -L$(BOOST_PYTHON_LIB) -lboost_python37 -lboost_numpy37

.PHONY : export_LD_LIBRARY_PATH print_LD_LIBRATY_PATH make_ldd_check clean

#$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $<
wrap: wrap.o
	$(CXX) -L/home/burmist/home2/training/032_boost/boost_1_76_0_install/lib -L/home/burmist/miniconda3/envs/boost_python/lib -L/home/burmist/miniconda3/envs/boost_python/lib/python3.7/config -Wl,-R -Wl,/home/burmist/home2/training/032_boost/boost_1_76_0_install/lib -Wl,-R -Wl,/home/burmist/miniconda3/envs/boost_python/lib -Wl,-R -Wl,/home/burmist/miniconda3/envs/boost_python/lib/python3.7/config -Wl,-rpath-link -Wl,/home/burmist/home2/training/032_boost/boost_1_76_0_install/lib -o $@ -Wl,--start-group $< -Wl,-Bstatic  -Wl,-Bdynamic -lboost_numpy37 -lboost_python37 -ldl -lpthread -lpython3.7 -lutil -Wl,--end-group -g

wrap.o: wrap.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

hello_ext.so: hello.o
	$(CXX) -o $@ -shared -Wl,--start-group hello.o $(LIB) -Wl,--end-group

log_gaussian_heap_ext.so: log_gaussian_heap.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -shared -Wl,--start-group log_gaussian_heap.o $(LIB) -Wl,--end-group

extending.so: extending.o
	$(CXX) -o $@ -shared -Wl,--start-group extending.o $(LIB) -Wl,--end-group

gaussian.so: gaussian.o
	$(CXX) -o $@ -shared -Wl,--start-group gaussian.o $(LIB) -Wl,--end-group

log_gaussian_boost_python.so: log_gaussian_boost_python.o
	$(CXX) -o $@ -shared -Wl,--start-group $< $(LIB) -Wl,--end-group

hello.o: hello.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

gaussian.o: gaussian.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

extending.o: extending.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

log_gaussian_heap.o: log_gaussian_heap.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

log_gaussian_boost_python.o: log_gaussian_boost_python.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

export_LD_LIBRARY_PATH:
	@echo export LD_LIBRARY_PATH=$(BOOST_PYTHON_LIB):$(LD_LIBRARY_PATH)

print_LD_LIBRATY_PATH:
	@echo LD_LIBRARY_PATH = $(value LD_LIBRARY_PATH)

make_ldd_check	:
	@echo 
	ldd hello_ext.so
	@echo 
	ldd extending.so

clean:
	rm -f *~
	rm -f .*~
	rm -f *.o
	rm -f *.so
	rm -f wrap
	rm -rf __pycache__
