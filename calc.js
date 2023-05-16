function AllCalculations(amplitude, numOfJ, Z=[], length=[], velocity=[], twoj=0, numOfIterations) {
    // let numOfIterations = 1000;
    let Vi = amplitude*1000;
    let Ii = Z[0]===0?Vi/(1e-15):Vi/Z[0]
    //for loop to create time to travel array
    let T = []
    for (let i = 0; i < numOfJ-1; i++) {
        T.push(length[i]/velocity[i])
    } 
    let tauF = [];
    let tauR = [];
    let rhoF = [];
    let rhoR = [];
    // initiate nested arrayes based on number of junctions
    let time = Array.from({ length: numOfJ }, () => [])
    let voltage = Array.from({ length: numOfJ }, () => [])
    let voltageL = Array.from({ length: numOfJ }, () => [])
    let voltageR = Array.from({ length: numOfJ }, () => [])
    let voltageTr = Array.from({ length: numOfJ }, () => [])
    let current = Array.from({ length: numOfJ }, () => [])
    let currentL = Array.from({ length: numOfJ }, () => [])
    let currentR = Array.from({ length: numOfJ }, () => [])
    let currentTr = Array.from({ length: numOfJ }, () => [])
    // for loop to calculate tau forward for every junction
    for (let i = 0; i < numOfJ; i++) {
        if(Z[i+1]===Infinity){
            tauF.push(2);
        }
        else if(Z[i+1].length){
            tauF.push(2/(Z[i]*((1/Z[i])+(1/Z[i+1][0])+((1/Z[i+1][1])||0))));
        }
        else{
            tauF.push((2*Z[i+1]/(Z[i]+Z[i+1])));
        }
        //get rho forward, tau and rho reverse for the same junction based on tau forward value with some relations between tau and rho
        rhoF.push(tauF[i]-1);
        tauR.push(1-rhoF[i]);
        rhoR.push(tauR[i]-1);
    }
    // get tau and rho for the current based on relations between current and voltage junction coeffecients
    let tauiF = tauR;
    let tauiR = tauF;
    let rhoiF = rhoR;
    let rhoiR = rhoF;
    let voltCoefficients = [tauF, tauR, rhoF,rhoR]
    let currentCoeffecients = [tauiF, tauiR, rhoiF, rhoiR]
    //----------------------middle volt and current----------------------------
    let tInitial = 0
    voltageR[0].push(Vi*tauF[0])
    currentR[0].push(Ii*tauiF[0])
    let m = 1
    while (m<=numOfJ-2) {
        tInitial = tInitial + T[m-1]
        time[m].push(tInitial)
        voltageL[m].push(voltageR[m-1][0]*rhoF[m])
        voltageR[m].push(voltageR[m-1][0]*tauF[m])
        voltage[m].push(voltageR[m][0])
        voltageTr[m].push(voltageR[m][0])
        //---------------------------
        currentL[m].push(currentR[m-1][0]*rhoiF[m])
        currentR[m].push(currentR[m-1][0]*tauiF[m])
        current[m].push(currentR[m][0])
        currentTr[m].push(currentR[m][0])
        let i = 1
        while (i<=numOfIterations) {
           if(m===1){
                time[m].push(time[m][i-1]+2*T[m-1]);

                voltageL[m].push((voltageL[m][i-1])*rhoR[m-1]*rhoF[m]);
                voltageR[m].push((voltageL[m][i-1])*rhoR[m-1]*tauF[m]);
                //voltage[m].push((voltage[m][voltage[m].length-1])+voltageR[m][voltageR[m].length-1]);
                voltageTr[m].push(voltageR[m][voltageR[m].length-1]);
                //----------------
                currentL[m].push((currentL[m][i-1])*rhoiR[m-1]*rhoiF[m]);
                currentR[m].push((currentL[m][i-1])*rhoiR[m-1]*tauiF[m]);
                //current[m].push((current[m][current[m].length-1])+currentR[m][currentR[m].length-1]);
                currentTr[m].push(currentR[m][currentR[m].length-1]);
           }else{
                time[m].push(time[m-1][i]+T[m-1]);

                voltageL[m].push((voltageR[m-1][i])*rhoF[m]);
                voltageR[m].push((voltageR[m-1][i])*tauF[m]);
                voltageTr[m].push(voltageR[m][voltageR[m].length-1]);
                //--------------------------------------
                currentL[m].push((currentR[m-1][i])*rhoiF[m]);
                currentR[m].push((currentR[m-1][i])*tauiF[m]);
                currentTr[m].push(currentR[m][currentR[m].length-1]);
           }

            let Tadd = T[m]
            let TauMultiply = 1
            //--------------------------
            let TauiMultiply = 1
            let j = 1
            while (j<= numOfJ-m-1) {
                if(j>1){
                    Tadd = Tadd +T[m-1+j]
                    TauMultiply = TauMultiply*tauF[j]*tauR[j]
                    Tadd = Tadd +T[m-1+j]
                    TauiMultiply = TauiMultiply*tauiF[j]*tauiR[j]
                
                } 
                time[m].push(time[m][i-1]+2*Tadd)
                voltageL[m].push(voltageR[m][i-1]*tauR[m]*rhoF[m+j]*TauMultiply)
                voltageR[m].push(voltageR[m][i-1]*rhoR[m]*rhoF[m+j]*TauMultiply)
                //voltage[m].push((voltage[m][voltage[m].length-1])+voltageL[m][voltageL[m].length-1])
                voltageTr[m].push(voltageL[m][voltageL[m].length-1])
                //-------------current
                currentL[m].push(currentR[m][i-1]*tauiR[m]*rhoiF[m+j]*TauiMultiply)
                currentR[m].push(currentR[m][i-1]*rhoiR[m]*rhoiF[m+j]*TauiMultiply)
                // current[m].push((current[m][current[m].length-1])+currentL[m][currentL[m].length-1])
                currentTr[m].push(currentL[m][currentL[m].length-1])
                j++
            }
            i++
        }
        m++
    }
    //-------------------------edge volt and current---------------------------
    let i = 1
    time[0].push(0)
    time[numOfJ-1].push(0)
    voltage[0].push(Vi*tauF[0])
    voltage[numOfJ-1].push(0)
    voltageTr[0].push(voltage[0][0])
    voltageTr[numOfJ-1].push(0)
    //------------------------------
    current[0].push(Ii*tauiF[0])
    current[numOfJ-1].push(0)
    currentTr[0].push(current[0][0])
    currentTr[numOfJ-1].push(0)
  //  let j = 2;
    let t01 = T[0]
    let t02 = T[numOfJ-2]
    while(i<=numOfIterations){
        // if(twoj===1){
        //     j = i+1;
        //     t01 = -T[0]
        //     t02 = -T[numOfJ-2]
        // } else {
        //     j=i;
        // }
            time[0].push(time[1][i-1]+t01)
            time[numOfJ-1].push(time[numOfJ-2][i-1]+t02)
            voltageTr[0].push((voltageL[1][i-1])*tauR[0])
            voltageTr[numOfJ-1].push((voltageR[numOfJ-2][i-1])*tauF[numOfJ-1])
            currentTr[0].push((currentL[1][i-1])*tauiR[0])
            currentTr[numOfJ-1].push((currentR[numOfJ-2][i-1])*tauiF[numOfJ-1])
            i++
        }

//---------------------------remove repeated values---------------------
    let newVoltage = Array.from({ length: numOfJ }, () => [])
    let newCurrent = Array.from({ length: numOfJ }, () => [])
    let newTime = Array.from({ length: numOfJ }, () => [])
        let temp;
        
        function converter(arrV,arrC,arrT) {
            for(let i=0; i<arrT.length; i++) {
          
              for (let j=i+1; j<arrT.length; j++) {
          
                if(arrT[i] > arrT[j]) {
          
                  temp = arrT[i]
                  arrT[i] = arrT[j]
                  arrT[j] = temp
                  temp = arrV[i]
                  arrV[i] = arrV[j]
                  arrV[j] = temp
                  temp = arrC[i]
                  arrC[i] = arrC[j]
                  arrC[j] = temp
                }
              }
            }
            return [arrV, arrC, arrT]
          }
        time.map((arr, index)=>converter(voltageTr[index], currentTr[index], time[index]))
        //---------------------accumiltate V and I transmitted--------------------
        for (let i = 0; i < numOfJ; i++) {
            for (let j = 1; j < time[i].length; j++) {
                voltage[i].push(voltage[i][j-1]+voltageTr[i][j])
                current[i].push(current[i][j-1]+currentTr[i][j])
            }
        }
        //-----------------get final value if repeated-------------------
        for (let i = 0; i < numOfJ; i++) {
            for (let j = time[i].length - 1; j >= 0; j--) {
                    if(isNaN(voltage[i][j])){
                        continue;
                    }else if(!newTime[i].includes(+time[i][j].toFixed(5))){
                        newTime[i].unshift(+time[i][j].toFixed(5))
                        newVoltage[i].unshift(+voltage[i][j].toFixed(5))
                        newCurrent[i].unshift(+current[i][j].toFixed(5))
                        }
                    }
            }

        if(twoj){
            return [
                    [newVoltage[0],newVoltage[2]],
                    [newCurrent[0],newCurrent[2]], 
                    [newTime[0],newTime[2]], 
                    [[tauF[0],tauF[2]], [tauR[0],tauR[2]], [rhoF[0],rhoF[2]],[rhoR[0],rhoR[2]]], 
                    [[tauiF[0],tauiF[2]], [tauiR[0],tauiR[2]], [rhoiF[0],rhoiF[2]],[rhoiR[0],rhoiR[2]]]
                ]
        }
        return [newVoltage, newCurrent, newTime, voltCoefficients, currentCoeffecients]
        // console.dir(time, {'maxArrayLength': null});
}
// console.log(AllCalculations(.5, 3, [200, 400, 40, Infinity], [450,300], [300,150]));

export {
AllCalculations
}