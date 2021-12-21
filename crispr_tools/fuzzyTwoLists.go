package main

/*
Writes pairs of sequences that differ by N characters, if there's
exactly ONE such sequence in the library file. Ignores exact matches
*/

import (
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)



func isClose(qSeq, libSeq string, reqMM int) bool {
	// from https://github.com/agnivade/levenshtein/blob/master/levenshtein.go
	// I'm making assumptions about the strings going in and removing safety checks

	if len(qSeq) != len(libSeq) {
		return false
	}
	mismatches := 0
	for i, _ := range qSeq {
		if qSeq[i] != libSeq[i] {
			mismatches++
		}
		if mismatches > reqMM {
			return false
		}
	}
	if mismatches == 1 {
		return true
	}
	return false
}


func check(e error) {
	if e != nil {
		panic(e)

	}
}

func loadStringArray(path string) []string {
	libFile, err := ioutil.ReadFile(path)
	check(err)
	libstr := string(libFile)
	lib := strings.Split(libstr, "\n")

	return lib
}

/* Compare a seq from file against every other in another.
Return the indicies of the sufficiently close
*/
func oneAgainstAll(query string, library *[]string, reqMM int) []int {
	var matches []int
	for i, libSeq := range *library {
		if isClose(query, libSeq, reqMM) {
			matches = append(matches, i)
		}
	}
	return matches
}


/* write pairs of sequences where the query sequence is  */
func compareFiles(queryFile string, libraryFile string, outFilePath string, reqMM int) {
	qArray := loadStringArray(queryFile)
	lArray := loadStringArray(libraryFile)

	// check we can write the out file before we do anything else
	fff, err := os.Create(outFilePath)
	check(err)
	err2 := fff.Close()
	check(err2)

	closeEnough := make(map[string]string, len(qArray))
	var waitGroup sync.WaitGroup

	semaphore := make(chan int, 1) //for mutex
	for _, query := range qArray {
		waitGroup.Add(1)
		go func(query string) {
			defer waitGroup.Done()
			matches := oneAgainstAll(query, &lArray, reqMM)

			if len(matches) == 1 {
				match := lArray[matches[0]]
				// lock, mutate, release
				semaphore <- 1
				closeEnough[query] = match
				<-semaphore
			}


		}(query)
	}
	waitGroup.Wait()
	outFile, err := os.Create(outFilePath)
	check(err)
	defer outFile.Close()

	var writeErr error
	nMatched := 0
	for wrongSeq, rightSeq := range closeEnough {
		nMatched++
		_, writeErr = outFile.WriteString(wrongSeq+"\t"+rightSeq+"\n")
	}
	check(writeErr)
	err4 := outFile.Close()
	check(err4)
	fmt.Println("Query seqs matched to similar:", nMatched, "of", len(qArray))
}

func main() {
	fmt.Println("fuzzyTwoLists.go Build 1003")

	// Use this when debugging from inside GoLand
	//os.Args = os.Args[1:len(os.Args)]
	if len(os.Args) < 4 {
		fmt.Println("fuzzyTwoLists QUERY_FILE LIBRARY_FILE OUT_FILE [DISTANCE default:1] [MAXPROCS]")
		os.Exit(1)
	}

	start := time.Now()
	qf := os.Args[1]
	lf := os.Args[2]
	outf := os.Args[3]

	var maxDist = 1
	var err error
	if len(os.Args) > 4 {
		maxDist, err = strconv.Atoi(os.Args[4])
		check(err)
	}

	// limit the number of processors
	if len(os.Args) > 5 {
		procs, err := strconv.Atoi(os.Args[5])
		check(err)
		maxProc := runtime.GOMAXPROCS(0)
		if procs < maxProc {
			runtime.GOMAXPROCS(procs)
		}
	}

	compareFiles(qf, lf, outf, maxDist)
	elapsed := time.Since(start)
	fmt.Println("Seconds elapsed:", int(elapsed.Seconds()))

}