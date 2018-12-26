TARGET = unit

CXXFLAGS = -Wall -std=c++14 -g

SRCS = $(TARGET).cpp
OBJS = $(SRCS:.cpp=.o)
HEADS = $(TARGET).h

$(TARGET): $(OBJS) $(HEADS)
	$(CXX) -o $(TARGET) $(OBJS) -lgmp -lgmpxx

all: $(TARGET)

run: all
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

