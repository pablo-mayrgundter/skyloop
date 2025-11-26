import { useFrame } from '@react-three/fiber'
import { useMemo, useRef } from 'react'
import { Mesh, Vector3 } from 'three'
import { LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle, LOOP_RADIUS_KM } from '../../config/simulation'
import { useSimulationStore } from '../../state/simulationStore'

const LoopGliders = () => {
  const gliderRef = useRef<Mesh>(null)
  const time = useSimulationStore((state) => state.time)
  const simSpeed = useSimulationStore((state) => state.simSpeed)

  const path = useMemo(() => {
    const pts: Vector3[] = []
    const segments = 64
    for (let i = 0; i <= segments; i += 1) {
      pts.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, i / segments, LOOP_RADIUS_KM))
    }
    return pts
  }, [])

  useFrame(() => {
    if (!gliderRef.current) return
    const loopTime = (time * simSpeed) % path.length
    const idx = Math.floor(loopTime)
    const t = loopTime - idx
    const current = path[idx % path.length]
    const next = path[(idx + 1) % path.length]
    const position = current.clone().lerp(next, t)
    gliderRef.current.position.copy(position)
    gliderRef.current.lookAt(next)
  })

  // Glider: 0.001 km = 1 meter radius, 0.005 km = 5 meters length (human-scale capsule)
  return (
    <mesh ref={gliderRef}>
      <capsuleGeometry args={[0.001, 0.005, 8, 16]} />
      <meshStandardMaterial color="#ffdf6b" emissive="#ffbe3b" emissiveIntensity={0.8} />
    </mesh>
  )
}

export default LoopGliders
