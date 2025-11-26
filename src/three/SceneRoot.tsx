import {Environment} from '@react-three/drei'
import {useMemo} from 'react'
import {Color, DoubleSide, MeshStandardMaterial} from 'three'
import {useSimulationStore} from '../state/simulationStore'
import BirdsEyeRig from './cameras/BirdsEyeRig'
import FirstPersonRig from './cameras/FirstPersonRig'
import Globe from './Globe'
import LoopGliders from './loop/LoopGliders'
import LoopSections from './loop/LoopSections'
import LoopTethers from './loop/LoopTethers'
import LoopTube from './loop/LoopTube'

const Atmosphere = () => {
  const material = useMemo(() => new MeshStandardMaterial({
    color: new Color('#4d74ff'),
    transparent: true,
    opacity: 0.08,
    side: DoubleSide
  }), [])
  return (
    <mesh>
      <sphereGeometry args={[0.07, 48, 48]} />
      <primitive object={material} />
    </mesh>
  )
}

const SceneRoot = () => {
  const mode = useSimulationStore((s) => s.mode)
  const {showAtmosphere, showGliders, showSections, showTethers} = useSimulationStore((s) => ({
    showAtmosphere: s.showAtmosphere,
    showGliders: s.showGliders,
    showSections: s.showSections,
    showTethers: s.showTethers
  }))
  return (
    <>
      {mode === 'firstPerson' ? <FirstPersonRig /> : <BirdsEyeRig />}
      <Globe />
      <LoopTube />
      {showSections && <LoopSections />}
      {showTethers && <LoopTethers />}
      {showGliders && <LoopGliders />}
      {showAtmosphere && <Atmosphere />}
      <Environment preset="sunset" />
    </>
  )
}

export default SceneRoot
