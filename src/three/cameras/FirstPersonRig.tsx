import { PerspectiveCamera } from '@react-three/drei'
import { useFrame } from '@react-three/fiber'
import { useMemo, useRef } from 'react'
import { PerspectiveCamera as PerspectiveCameraImpl, Vector3 } from 'three'
import { LOS_ANGELES, SAN_FRANCISCO } from '../../config/simulation'
import { useSimulationStore } from '../../state/simulationStore'
import { buildGreatCircleArc } from '../Globe'

const FirstPersonRig = () => {
  const cameraRef = useRef<PerspectiveCameraImpl>(null)
  const path = useMemo(() => buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, 512), [])
  const simSpeed = useSimulationStore((s) => s.simSpeed)
  const time = useSimulationStore((s) => s.time)

  useFrame((state, delta) => {
    const progress = ((time + delta * simSpeed) * 0.02) % 1
    const idx = Math.floor(progress * (path.length - 1))
    const next = path[idx]
    const lookAhead = path[(idx + 2) % path.length] ?? next.clone().add(new Vector3(0.01, 0, 0))
    if (cameraRef.current && next) {
      cameraRef.current.position.copy(next)
      cameraRef.current.lookAt(lookAhead)
      state.camera.updateProjectionMatrix()
    }
  })

  return <PerspectiveCamera ref={cameraRef} fov={60} near={0.001} far={100} makeDefault position={[0, 0.01, 0]} />
}

export default FirstPersonRig
