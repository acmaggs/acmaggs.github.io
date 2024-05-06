comm -3 <(find . | sort) <(cd ../crawl/acmaggs.github.io/ ; find . | sort) |grep -v git
