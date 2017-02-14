CXXFLAGS = -std=c++17 -O3 -march=native -g
CXX = g++

SRC = src/
BUILD = build/
COMPONENT_DIR = src/components/
COMPONENT_OUTPUT_DIR = build/components/

maxclique: $(SRC)/main.cpp $(COMPONENT_OUTPUT_DIR)/libhpx_incumbent_component.so $(COMPONENT_OUTPUT_DIR)/libhpx_workqueue_component.so $(BUILD)/DimacsParser.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(BUILD)/DimacsParser.o `pkg-config --cflags --libs hpx_application` -I$(COMPONENT_DIR) -DHPX_APPLICATION_NAME=incumbentTest -L$(COMPONENT_OUTPUT_DIR) -lhpx_incumbent_component -lhpx_iostreams -lhpx_workqueue_component

$(BUILD)/%.o: $(SRC)/%.cpp $(SRC)/%.hpp
	mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(COMPONENT_OUTPUT_DIR)/libhpx_incumbent_component.so: $(COMPONENT_DIR)/incumbent_component.cpp $(COMPONENT_DIR)/incumbent_component.hpp
	mkdir -p $(COMPONENT_OUTPUT_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< `pkg-config --cflags --libs hpx_component` -DHPX_COMPONENT_NAME=incumbent_component

$(COMPONENT_OUTPUT_DIR)/libhpx_workqueue_component.so: $(COMPONENT_DIR)/workqueue_component.cpp $(COMPONENT_DIR)/workqueue_component.hpp
	mkdir -p $(COMPONENT_OUTPUT_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $< `pkg-config --cflags --libs hpx_component` -DHPX_COMPONENT_NAME=workqueue_component
