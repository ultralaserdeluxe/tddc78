.PHONY: all
all: lab1 lab2 lab3 lab5

.PHONY: lab1
lab1:
	cd lab1 && $(MAKE)

.PHONY: lab2
lab2:
	cd lab2 && $(MAKE)

.PHONY: lab3
lab3:
	cd lab3 && $(MAKE)

.PHONY: lab5
lab5:
	cd lab5 && $(MAKE)

.PHONY: clean
clean:
	cd lab1 && $(MAKE) clean
	cd lab2 && $(MAKE) clean
	cd lab3 && $(MAKE) clean
	cd lab5 && $(MAKE) clean
