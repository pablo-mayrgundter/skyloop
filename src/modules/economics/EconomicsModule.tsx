const EconomicsModule = () => {
  return (
    <div className="controls-grid">
      <h3>Economics</h3>
      <p>Placeholder economic model exploring ticket pricing and opex/capex.</p>
      <ul>
        <li>Ticket price target: $120</li>
        <li>Daily passengers: 40,000</li>
        <li>Payback: ~9 years (rough cut)</li>
      </ul>
      <p>TODO: incorporate Monte Carlo uncertainty and demand modeling.</p>
    </div>
  )
}

export default EconomicsModule
