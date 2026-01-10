// A simple CLI for testing the LevelDB implementation.
package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	leveldb "leveldb-go"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: leveldb-cli <dbpath>")
		os.Exit(1)
	}

	dbPath := os.Args[1]
	db, err := leveldb.Open(dbPath, nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to open database: %v\n", err)
		os.Exit(1)
	}
	defer db.Close()

	fmt.Println("LevelDB CLI - Commands: put <key> <value>, get <key>, delete <key>, scan, quit")

	scanner := bufio.NewScanner(os.Stdin)
	fmt.Print("> ")
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			fmt.Print("> ")
			continue
		}

		parts := strings.SplitN(line, " ", 3)
		cmd := strings.ToLower(parts[0])

		switch cmd {
		case "put":
			if len(parts) < 3 {
				fmt.Println("Usage: put <key> <value>")
			} else {
				key := parts[1]
				value := parts[2]
				if err := db.Put([]byte(key), []byte(value)); err != nil {
					fmt.Fprintf(os.Stderr, "Error: %v\n", err)
				} else {
					fmt.Println("OK")
				}
			}

		case "get":
			if len(parts) < 2 {
				fmt.Println("Usage: get <key>")
			} else {
				key := parts[1]
				value, err := db.Get([]byte(key))
				if err == leveldb.ErrNotFound {
					fmt.Println("(not found)")
				} else if err != nil {
					fmt.Fprintf(os.Stderr, "Error: %v\n", err)
				} else {
					fmt.Printf("%s\n", value)
				}
			}

		case "delete", "del":
			if len(parts) < 2 {
				fmt.Println("Usage: delete <key>")
			} else {
				key := parts[1]
				if err := db.Delete([]byte(key)); err != nil {
					fmt.Fprintf(os.Stderr, "Error: %v\n", err)
				} else {
					fmt.Println("OK")
				}
			}

		case "scan":
			it := db.NewIterator()
			it.SeekToFirst()
			count := 0
			for it.Valid() {
				fmt.Printf("  %s -> %s\n", it.Key(), it.Value())
				count++
				it.Next()
			}
			it.Close()
			fmt.Printf("(%d entries)\n", count)

		case "quit", "exit", "q":
			fmt.Println("Goodbye!")
			return

		case "help":
			fmt.Println("Commands:")
			fmt.Println("  put <key> <value>  - Store a key-value pair")
			fmt.Println("  get <key>          - Retrieve a value")
			fmt.Println("  delete <key>       - Delete a key")
			fmt.Println("  scan               - Scan all entries")
			fmt.Println("  quit               - Exit the CLI")

		default:
			fmt.Printf("Unknown command: %s (type 'help' for commands)\n", cmd)
		}

		fmt.Print("> ")
	}
}

