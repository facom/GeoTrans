BRANCH=$(shell bash .getbranch)

branch:
	@echo "You are in: $(PWD)"
	@echo $(BRANCH)

clean:
	find . -name "*~" -exec rm -rf {} \;
	find . -name "#*#" -exec rm -rf {} \;
	find . -name "*.pyc" -exec rm -rf {} \;

cleandata:
	rm -rf confs/* tmp/* plots/*

cleanall:clean cleandata

commit:
	@echo "Commiting..."
	@-git commit -am "Commit"
	@-git push origin $(BRANCH)

pull:
	@echo "Pulling latest version..."
	@-git reset --hard HEAD
	@-git pull origin $(BRANCH)
