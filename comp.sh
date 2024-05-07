#compare two trees, the first crawled from ../crawl by wget, the second the local files, to find files that are not used
comm -3 <(find . | sort) <(cd ../crawl/acmaggs.github.io/ ; find . | sort) |grep -v git|grep -v DS_Store
