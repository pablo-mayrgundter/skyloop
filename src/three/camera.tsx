import { useFrame, useThree } from '@react-three/fiber'
import { useRef } from 'react'
import { PerspectiveCamera } from 'three'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO } from '../config/simulation'
import { useSimulationStore } from '../state/simulationStore'

export const useFirstPersonCamera = () => {
  const camera = useThree((state) => state.camera as PerspectiveCamera)
  const { time, simSpeed } = useSimulationStore((state) => ({
    time: state.time,
    simSpeed: state.simSpeed
  }))
  useFrame(() => {
    const loopTime = (time * simSpeed) % 300
    const t = (loopTime % 60) / 60
    const position = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t)
    camera.position.copy(position.clone().multiplyScalar(1.02))
    camera.lookAt(0, 0, 0)
  })
}

export const FirstPersonRig = () => {
  useFirstPersonCamera()
  return null
}

export const BirdsEyeRig = ({ overlays }: { overlays: Record<string, boolean> }) => {
  const camera = useThree((state) => state.camera as PerspectiveCamera)
  const targetRef = useRef({ x: 0, y: 0, z: 0 })

  useFrame(() => {
    targetRef.current = { x: 0, y: 0, z: 0 }
    camera.position.lerp({ x: 0, y: 3.4, z: 2.8 }, 0.05)
    camera.lookAt(targetRef.current.x, targetRef.current.y, targetRef.current.z)
  })

  return (
    <group>
      {overlays.showWind && (
        <mesh position={[0, 1.6, 0]}>
          <coneGeometry args={[0.05, 0.2, 6]} />
          <meshStandardMaterial color="#74f0ed" transparent opacity={0.6} />
        </mesh>
      )}
    </group>
  )
}
