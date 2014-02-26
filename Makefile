SRCS = main.cpp algorithm.cpp class.cpp common.cpp const.cpp 
OBJS = $(SRCS:.c=.cpp)
EXECUTABLE = NURD

all:
	g++ -O3 -o $(EXECUTABLE) $(SRCS)

.PHONY: clean

clean:
	@echo "cleaning project"
	-rm $(OBJ_DIR)/*.o
	@echo "clean completed"
