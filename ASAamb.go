package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/align/matrix"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

func getFasta(fn string) (seq.Sequence, error) {
	fasta_file, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	defer fasta_file.Close()
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.Protein)
	reader := fasta.NewReader(fasta_file, t)
	seq, _ := reader.Read()
	return seq, nil
}

func getSurf(fn string) ([]float64, error) {
	var s []float64
	var f float64
	surf_file, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	defer surf_file.Close()
	scanner := bufio.NewScanner(surf_file)
	scanner.Scan()

	scanner.Split(bufio.ScanWords)
	scanner.Scan()
	for scanner.Scan() {
		f, _ = strconv.ParseFloat(scanner.Text(), 64)
		if f < 0.0 {
			f = 0.0
		}
		s = append(s, f)
	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, "reading standard input:", err)
	}
	// fmt.Println(s)

	return s, err
}

// func getAttr(fn string) (seq.Sequence, seq.Sequence, error) {
// 	fasta_file, err := os.Open(fn)
// 	if err != nil {
// 		fmt.Println("Erro ao ler o arquivo", err)
// 	}
// 	var s []alphabet.Letter
// 	t := linear.NewSeq("", s, alphabet.Protein)
// 	reader := fasta.NewReader(fasta_file, t)
// 	seq, _ := reader.Read()
// 	// if reader.Next() != nil {
// 	// 	fmt.Println("PAU NO NEXT")
// 	// }
// 	attr, _ := reader.Read()
// 	//fmt.Println("Read -> ", seq.Alphabet())
// 	return seq, attr, nil
// }

func alinha(seq_aa, seqpdb_aa seq.Sequence, surf []float64) (string, string) {
	nw := align.NWAffine{
		Matrix:  matrix.MATCH,
		GapOpen: -1,
	}

	aln, err := nw.Align(seq_aa, seqpdb_aa)
	var f_aa [2]alphabet.Slice
	if err == nil {
		f_aa = align.Format(seq_aa, seqpdb_aa, aln, '_')
		// f_ss = align.Format(seq_aa, ss_ss, aln, '_')
	} else {
		fmt.Println("O ERRO E:", err)
	}
	return fmt.Sprint(f_aa[0]), fmt.Sprint(f_aa[1])
}

func main() {

	fastaFName := flag.String("fa", "", "fasta file")
	fastapdbFName := flag.String("fapdb", "", "fasta pdb file")
	surfFName := flag.String("surf", "", "surf file")

	flag.Parse()

	var seq_aa, seqpdb_aa seq.Sequence
	var surf []float64
	var err error

	// Le a sequencia fasta
	if *fastaFName == "" {
		fmt.Println("É necessario o arquivo fasta")
		return
	} else {
		// seq_aa, err = getFasta("/home/jgcarvalho/sync/ZZpred/pisces/fasta/1A1XA.fa")
		seq_aa, err = getFasta(*fastaFName)
		if err != nil {
			fmt.Println("Erro no processamento do Fasta", err)
		}
	}

	if *fastapdbFName == "" {
		fmt.Println("É necessario o arquivo Fasta do PDB")
		return
	} else {
		// dssp_aa, dssp_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/dssp_chain_fasta/1A1XA_dssp.fasta")
		seqpdb_aa, err = getFasta(*fastapdbFName)
		if err != nil {
			fmt.Println("Erro no processamento do Fasta PDB", err)
		}
	}

	if *surfFName == "" {
		fmt.Println("É necessario o arquivo de ASA")
		return
	} else {
		// stride_aa, stride_ss, err = getAttr("/home/jgcarvalho/sync/ZZpred/pisces/stride_chain_fasta/1A1XA_stride.fasta")
		surf, err = getSurf(*surfFName)
		if err != nil {
			fmt.Println("Erro no processamento do arquivo de ASA", err)
		}
	}

	// surf := []float64{0.0, 0.0}
	// s1, s2 := alinha(seq_aa, seqpdb_aa, surf)
	_, s2 := alinha(seq_aa, seqpdb_aa, surf)

	// fmt.Println(s1)
	// fmt.Println(s2)
	pep3 := strings.Split(s2, "")
	// fmt.Println(pep3)
	tmp := make([]float64, len(pep3), len(pep3))
	count := 0
	for i, v := range pep3 {
		if v == "_" {
			tmp[i] = 0.0
		} else {
			tmp[i] = surf[count]
			count++
		}
	}
	surf = tmp
	// fmt.Println(surf)

	for i := 1; i < (len(pep3) - 1); i++ {
		if (pep3[i-1] != "_") && (pep3[i] != "_") && (pep3[i+1] != "_") {
			fmt.Printf("%s%s%s %f\n", pep3[i-1], pep3[i], pep3[i+1], surf[i])
		}
	}

	// Monta os alinhamentos
	// aa, dsspaa, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	// aa, _, dsspss := alinha(seq_aa, dssp_aa, dssp_ss)
	// // _, strideaa, stridess := alinha(seq_aa, stride_aa, stride_ss)
	// _, _, stridess := alinha(seq_aa, stride_aa, stride_ss)
	// // _, prossaa, prossss := alinha(seq_aa, pross_aa, pross_ss)
	// _, _, prossss := alinha(seq_aa, pross_aa, pross_ss)
	// // _, kaksiaa, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)
	// _, _, kaksiss := alinha(seq_aa, kaksi_aa, kaksi_ss)

	// fmt.Printf("PDBSEQUENCE %s\n", aa)
	// // fmt.Printf("DSSPSEQ     %s\n", dsspaa)
	// fmt.Printf("DSSPSECSTR  %s\n", dsspss)
	// // fmt.Printf("STRDSEQ     %s\n", strideaa)
	// fmt.Printf("STRDSECSTR  %s\n", stridess)
	// // fmt.Printf("PROSSEQ     %s\n", prossaa)
	// fmt.Printf("PROSSECSTR  %s\n", prossss)
	// // fmt.Printf("KAKSSEQ     %s\n", kaksiaa)
	// fmt.Printf("KAKSSECSTR  %s\n", kaksiss)

	// Escreve um arquivo com os dados de ss e gaps em um arquivo JSON
	// data := Data{
	// 	Seq:     aa,
	// 	Dssp:    dsspss,
	// 	Stride:  stridess,
	// 	Kaksi:   kaksiss,
	// 	Pross:   prossss,
	// 	Dssp3:   dssp3(dsspss),
	// 	Stride3: stride3(stridess),
	// 	Kaksi3:  kaksi3(kaksiss),
	// 	Pross3:  pross3(prossss),
	// }
	// b, err := json.MarshalIndent(data, "", "  ")
	// if err != nil {
	// 	fmt.Println("error:", err)
	// }
	// fmt.Println(string(b))
}
