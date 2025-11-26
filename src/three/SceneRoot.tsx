import { Suspense } from 'react'
import { Environment, Grid, Html } from '@react-three/drei'
import Globe from './Globe'
import LoopSections from './loop/LoopSections'
import LoopTethers from './loop/LoopTethers'
import LoopTube from './loop/LoopTube'
import LoopGliders from './loop/LoopGliders'
import { BirdsEyeRig, FirstPersonRig } from './camera'
import { useSimulationStore } from '../state/simulationStore'

const SceneRoot = () => {
  const mode = useSimulationStore((state) => state.mode)
  const overlays = useSimulationStore((state) => ({
    showAtmosphere: state.showAtmosphere,
    showWind: state.showWind,
    showTethers: state.showTethers,
    showSections: state.showSections,
    showGliders: state.showGliders,
    showDayNight: state.showDayNight
  }))

  return (
    <Suspense fallback={<Html center>Loading assets...</Html>}>
      <Environment preset="sunset" />
      <Grid infiniteGrid cellSize={0.1} sectionSize={2} sectionColor="#1a2a42" />
      <Globe overlays={overlays} />
      {overlays.showSections && <LoopSections />}
      {overlays.showTethers && <LoopTethers />}
      <LoopTube />
      {overlays.showGliders && <LoopGliders />}
      {mode === 'firstPerson' ? <FirstPersonRig /> : <BirdsEyeRig overlays={overlays} />}
    </Suspense>
  )
}

export default SceneRoot
