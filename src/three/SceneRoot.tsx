import { Suspense, useMemo } from 'react'
import { Environment, Html } from '@react-three/drei'
import Globe from './Globe'
import Grid from './Grid'
import LoopSections from './loop/LoopSections'
import LoopTethers from './loop/LoopTethers'
import LoopTube from './loop/LoopTube'
import LoopGliders from './loop/LoopGliders'
import { BirdsEyeRig, FirstPersonRig, FollowRig } from './camera'
import { useSimulationStore } from '../state/simulationStore'

const SceneRoot = () => {
  const mode = useSimulationStore((state) => state.mode)
  const showAtmosphere = useSimulationStore((state) => state.showAtmosphere)
  const showWind = useSimulationStore((state) => state.showWind)
  const showTethers = useSimulationStore((state) => state.showTethers)
  const showSections = useSimulationStore((state) => state.showSections)
  const showGliders = useSimulationStore((state) => state.showGliders)
  const showDayNight = useSimulationStore((state) => state.showDayNight)
  
  const overlays = useMemo(
    () => ({
      showAtmosphere,
      showWind,
      showTethers,
      showSections,
      showGliders,
      showDayNight
    }),
    [showAtmosphere, showWind, showTethers, showSections, showGliders, showDayNight]
  )

  return (
    <Suspense fallback={<Html center>Loading assets...</Html>}>
      <Environment preset="sunset" />
      {/* Grid disabled for now - was creating plane effect when viewed from space */}
      {/* <Grid infiniteGrid cellSize={1000} sectionSize={10000} sectionColor="#1a2a42" /> */}
      <Globe overlays={overlays} />
      {overlays.showSections && <LoopSections />}
      {overlays.showTethers && <LoopTethers />}
      <LoopTube />
      {overlays.showGliders && <LoopGliders />}
      {mode === 'firstPerson' ? (
        <FirstPersonRig />
      ) : mode === 'follow' ? (
        <FollowRig />
      ) : (
        <BirdsEyeRig overlays={overlays} />
      )}
    </Suspense>
  )
}

export default SceneRoot
