package main

import (
	"C"
	"fmt"
	"math/rand"
	"net"
	"os"

	"bufio"
	"strconv"

	// "strconv"
	// "strings"
	"sync"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v3/dckks"
	"github.com/tuneinsight/lattigo/v3/drlwe"

	// "github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)
import (
	"strings"
	"time"

	//

	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

type SafeStringArray struct {
	mu          sync.Mutex
	stringArray []string
	numContents int
}

var publicKeyStr string
var cpkgSharesQStr []string
var cpkgSharesPStr []string
var cpkgWrite []sync.Mutex
var pkStrWrite []sync.Mutex
var done []sync.Mutex

// func (c *SafeStringArray) Update(str string) {
// 	c.mu.Lock()
// 	c.stringArray[c.numContents] = str
// 	c.numContents++
// 	c.mu.Unlock()
// }
func polyCoeffsEncode(coeffs [][]uint64) string {
	res := ""
	for dimZeroCounter := range coeffs {
		tmp := coeffs[dimZeroCounter]
		tmpRes := arrayToString(tmp)
		res += tmpRes + " "
	}
	return res
}

func arrayToString(arr []uint64) string {
	init := strconv.FormatUint(arr[0], 10)
	for cntr := 1; cntr < len(arr); cntr++ {
		init = init + "." + strconv.FormatUint(arr[cntr], 10)
	}
	return init
}
func polyCoeffsDecode(str string) [][]uint64 {
	//fmt.Println("I am called")
	polyCoeffsArr := strings.Split(str, " ")
	polyCoeffsArr = polyCoeffsArr[0 : len(polyCoeffsArr)-1]
	//fmt.Println("I am allocating")
	res := make([][]uint64, len(polyCoeffsArr))
	//fmt.Println("I allocated")
	for polyCounter := range res {
		tmp := stringToArray(polyCoeffsArr[polyCounter])
		res[polyCounter] = tmp
	}
	//fmt.Println("I am returning")
	return res
}
func stringToArray(strIn string) []uint64 {
	//fmt.Println("string to arr called")
	//fmt.Println(strIn[0:20])
	elements := strings.Split(strIn, ".")
	//fmt.Println("string to arr allocating")
	resLoc := make([]uint64, len(elements))
	//fmt.Println("string to arr allocated")
	for counter := range resLoc {
		resLoc[counter], _ = strconv.ParseUint(elements[counter], 10, 64)
	}
	//fmt.Println("string to arr returned")
	return resLoc
}

func handleClientPhase1(connections []net.Conn, idx int) {

	////// transmit the index of the client and number of peers
	fmt.Println("Handling client")
	numPeers := len(connections)
	conn := connections[idx]
	conn.Write([]byte(strconv.Itoa(idx) + " " + strconv.Itoa(numPeers) + "\n"))
	fmt.Println("Sending client index")

	/////// Receive shares and generate collective public key
	message, err := bufio.NewReader(conn).ReadString('\n')
	fmt.Println("Receiving client pk shares")
	if err != nil {
		panic(err)
	}
	sharesArr := strings.Split(message, "&")
	fmt.Println(len(sharesArr))
	sharesQArr := sharesArr[0]
	sharesPArr := sharesArr[1]

	cpkgSharesQStr[idx] = sharesQArr
	cpkgSharesPStr[idx] = sharesPArr
	cpkgWrite[idx].Unlock()
	pkStrWrite[idx].Lock()
	conn.Write([]byte(publicKeyStr))
	done[idx].Unlock()
}

func ServerSetup(serverAddress string, numPeers int) {

	l, err := net.Listen("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error listening:", err.Error())
		os.Exit(1)
	}
	// Close the listener when the application closes.
	defer l.Close()
	// The array to save the address of clients
	clientIPs := make([]net.Conn, numPeers)
	for cntr := 0; cntr < numPeers; cntr++ {
		// Listen for an incoming connection.
		c, err := l.Accept()
		if err != nil {
			fmt.Println("Error connecting:", err.Error())
			return
		}
		// Print client connection address.
		//fmt.Println("Client " + c.RemoteAddr().String() + " connected.")
		// Add connection to list
		clientIPs[cntr] = c
	}
	cpkgShares := make([]*drlwe.CKGShare, numPeers)
	cpkgSharesQStr = make([]string, numPeers)
	cpkgSharesPStr = make([]string, numPeers)
	cpkgWrite = make([]sync.Mutex, numPeers)
	done = make([]sync.Mutex, numPeers)
	pkStrWrite = make([]sync.Mutex, numPeers)
	for peerIdx := range cpkgWrite {
		pkStrWrite[peerIdx].Lock()
		cpkgWrite[peerIdx].Lock()
		done[peerIdx].Lock()

	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientPhase1(clientIPs, idx)
	}

	for peerIdx := range cpkgWrite {
		cpkgWrite[peerIdx].Lock()
	}

	//////////////////////////
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	params, err := ckks.NewParametersFromLiteral(ckksParams)
	ckg := dckks.NewCKGProtocol(params)
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	ckgCombined := ckg.AllocateShare()
	// ringQP, _ := ring.NewRing(params.N(), append(params.Q(), params.P()...))
	pk := ckks.NewPublicKey(params)
	var _ drlwe.CollectivePublicKeyGenerator = dckks.NewCKGProtocol(params)
	// cpkgShares := make([]dckks.CKGShare, numPeers)

	for peerIdx := range clientIPs {
		fmt.Println(clientIPs[peerIdx])
		coeffsQ := polyCoeffsDecode(cpkgSharesQStr[peerIdx])
		coeffsP := polyCoeffsDecode(cpkgSharesPStr[peerIdx])
		poly := ringqp.Poly{}
		poly.P.Coeffs = coeffsP
		poly.Q.Coeffs = coeffsQ
		// ringqp.Poly{{Q:coeffsQ,P:coeffsP}}
		// poly[0] = ring.Poly{Coeffs:coeffsQ}

		cpkgShares[peerIdx].Value = poly

		ckg.AggregateShare(cpkgShares[peerIdx], ckgCombined, ckgCombined)

	}
	crp := ckg.SampleCRP(prng)
	ckg.GenPublicKey(ckgCombined, crp, pk)

	publicKeyStr = ""
	// pkContent := pk.Value

	//////////////////////////
	// for itemIdx := range pkContent {
	publicKeyStr += polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "&" + polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "\n"
	// }
	// publicKeyStr += "\n"
	for peerIdx := range cpkgWrite {
		pkStrWrite[peerIdx].Unlock()
	}

	for peerIdx := range done {
		done[peerIdx].Lock()
		done[peerIdx].Unlock()
	}
}

func ClientSetup(serverAddress string) {

	//Build up communication
	// var ringPrime uint64 = 0x10000000001d0001
	// var ringPrimeP uint64 = 0xfffffffffffc001
	///////// Start the client and connect to the server.
	//startTime := time.Now()
	rand.Seed(time.Now().UTC().UnixNano())
	//fmt.Println("Connecting to", "tcp", "server", serverAddress)
	conn, err := net.Dial("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error connecting:", err.Error())
		os.Exit(1)
	}
	defer conn.Close()
	_, err = bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Connected")
	// moduli := &ckks.Moduli{Qi: []uint64{ringPrime}, Pi: []uint64{ringPrimeP}}
	// params, err := ckks.NewParametersFromLiteral(logDegree, moduli)
	// params.SetScale(scale)
	// params.SetLogSlots(logDegree - 1)
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	prng, _ := utils.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	params, err := ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	type Party struct {
		*dckks.CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}
	// minLevel, logBound, ok := dckks.GetMinimumLevelForBootstrapping(128, params.DefaultScale(), numPeers, params.Q())

	p := new(Party)
	p.CKGProtocol = dckks.NewCKGProtocol(params)
	p.s = ckks.NewSecretKey(params) //sk0Shards[0]
	fmt.Println("Gen sk")
	p.s1 = p.AllocateShare()

	crp := p.SampleCRP(prng)

	p.GenShare(p.s, crp, p.s1)
	fmt.Println("Gen pk")
	toSendString := ""
	fmt.Println("Gen string")
	toSendString += polyCoeffsEncode(p.s1.Value.Q.Coeffs)
	fmt.Println("Gen string 1")
	toSendString += "&"
	fmt.Println("Gen string 2")
	toSendString += polyCoeffsEncode(p.s1.Value.P.Coeffs)
	fmt.Println("Gen string 3")
	toSendString += "\n"
	fmt.Println("Gen sending strings")
	conn.Write([]byte(toSendString))
	fmt.Println("Sending pk share")

}

func Serverinference() {

}

func Clientinference() {

}

func main() {
	args := os.Args[1:]
	appType := args[0]
	fmt.Println(appType)
	// scale := float64(1 << 40)
	// robust := false
	// logDegree := uint64(13)
	numPeers := 2
	inLen := 20000
	if appType == "client" {
		inputs := make([]float64, inLen)
		for i := range inputs {
			inputs[i] = rand.Float64()
		}
		fmt.Println("Client setup")
		ClientSetup("localhost:8080")
		// clientPhase2(inputs, cpk, shamirShare, id, "localhost:8080", robust, logDegree, scale, 0.5)
	} else {
		fmt.Println("Server setup")
		ServerSetup("localhost:8080", numPeers)
		// serverPhase2("localhost:8080", numPeers, robust, 0.5, logDegree, scale, inLen)
	}
}
