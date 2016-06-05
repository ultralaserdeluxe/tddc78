.PHONY: all
all: blur threshold

.PHONY: blur
blur:
	cd blur	&& $(MAKE)

.PHONY: threshold
threshold:
	cd threshold && $(MAKE)

.PHONY: clean
clean:
	cd blur && $(MAKE) clean
	cd threshold && $(MAKE) clean
